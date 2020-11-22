# -*- coding: utf-8 -*-

# Copyright 2019 Zhenjun-Zhang
# email: zjzhang@ynao.ac.cn
# email: cnzjzhang@gmail.com
# Distributed under the terms of the YNAO license

import warnings
import ctypes
from ctypes import byref, POINTER, Structure
from ctypes import c_int, c_char, c_double, c_longdouble
from numpy.ctypeslib import ndpointer
from numpy import ndarray, array, zeros, asarray, asfarray

ll = ctypes.windll.LoadLibrary
_sofa = ll('H:/workF/sofa64.dll')


def has_function(funcname):
    """ Helper function that returns True if this particular release of |SOFA|
    provides the function *funcname*, and False otherwise. This is only the
    case with function names that are legal |SOFA| function names, wich means
    that calling ``has_function`` with a *funcname* that isn't known by any
    version of |SOFA| will raise :exc:`AttributeError`.

    """

    if not funcname in globals():
        raise AttributeError('%s does not know any function "named" %s' % \
                                (__package__, funcname))
    # convert 'funcname' to its equivalent SOFA name
    funcname = 'iau' + funcname[0].upper() + funcname[1:]
    return hasattr(_sofa, funcname)


#################################################
# SOFA Global Define
####################################################################
global DPI
DPI = 3.141592653589793238462643
global D2PI
D2PI = 6.283185307179586476925287
global DR2H  # Radians to hours
DR2H = 3.819718634205488058453210
global DR2S  # Radians to seconds
DR2S = 13750.98708313975701043156
global DR2D  # radians to degrees
DR2D = 57.29577951308232087679815
global DR2AS  # radians to acseconds
DR2AS = 206264.8062470963551564734
global DH2R  # Hours to radians
DH2R = 0.2617993877991494365385536
global DS2R  # Seconds to radinas
DS2R = 7.272205216643039903848712e-5
global DD2R  # Degrees to radians
DD2R = 1.745329251994329576923691e-2
global DAS2R # Arc seconds to radians
DAS2R = 4.848136811095359935899141e-6
#####################################################################################

#/* Milliarcseconds to radians */
DMAS2R =DAS2R / 1e3
#/* Arcseconds in a full circle */
TURNAS = 1296000.0
#/* Length of tropical year B1900 (days) */
DTY =365.242198781
#/* Seconds per day. */
DAYSEC =86400.0
#/* Days per Julian year */
DJY = 365.25
#/* Days per Julian century */
DJC = 36525.0
#/* Days per Julian millennium */
DJM = 365250.0
#/* Reference epoch (J2000.0), Julian Date */
DJ00 = 2451545.0
#/* Julian Date of Modified Julian Date zero */
DJM0 = 2400000.5
#/* Reference epoch (J2000.0), Modified Julian Date */
DJM00 = 51544.5
#/* 1977 Jan 1.0 as MJD */
DJM77 = 43144.0
#/* TT minus TAI (s) */
TTMTAI = 32.184
#/* Astronomical unit (m, IAU 2012) */
DAU = 149597870.7e3
#/* Speed of light (m/s) */
CMPS = 299792458.0
#/* Light time for 1 au (s) */
AULT = DAU/CMPS
#/* Speed of light (au per day) */
DC = DAYSEC/AULT
####################################################
#使用类来代替c语言中的结构体
class ASTROM(Structure):  
    _fields_ = [("pmt", c_double), #/* PM time interval (SSB, Julian years) */
                ("eb", c_double*3),# /* SSB to observer (vector, au) */
                ("eh", c_double*3),#/* Sun to observer (unit vector) */
                ("em", c_double),#/* distance from Sun to observer (au) */
                ("v", c_double*3),#/* barycentric observer velocity (vector, c) */
                ("bml", c_double),#/* sqrt(1-|v|^2): reciprocal of Lorenz factor */
                ("bpn",( c_double*3)*3),#/* bias-precession-nutation matrix */
                ("along", c_double),#/* longitude + s' + dERA(DUT) (radians) */
                ("phi", c_double),#/* geodetic latitude (radians) */
                ("xpl", c_double),#/* polar motion xp wrt local meridian (radians) */
                ("ypl", c_double),#/* polar motion yp wrt local meridian (radians) */
                ("sphi", c_double),#/* sine of geodetic latitude */
                ("cphi", c_double),#/* cosine of geodetic latitude */
                ("diurab", c_double),#/* magnitude of diurnal aberration vector */
                ("eral", c_double),#/* "local" Earth rotation angle (radians) */
                ("refa", c_double),#/* refraction constant A (radians) */
                ("refb", c_double)# /* refraction constant B (radians) */]
                ]  

class LDBODY(Structure):
    _fields_ = [("bm", c_double), #/* mass of the body (solar masses) */
                ("dl", c_double), #/* deflection limiter (radians^2/2) */
                ("pv",( c_double*3)*2) #/barycentric PV of the body (au, au/day) */
                ]  
####################################################################################
# SOFA Astrometry  TOOLs
#################################################################################



##############################################################################
'''
Star proper motion: update star catalog data for space motion, with
special handling to handle the zero parallax case
'''







################################################################
###########  ROUTINES Calendars
###############################################################
    
#iauCal2jd
try:
    _sofa.iauCal2jd.argtypes = [
    c_int,  # iy \year
    c_int,  # im \month
    c_int,  # idd \day
    POINTER(c_double),  # djm0\ MJD zero−point: always 2400000.5
    POINTER(c_double),  # djm \ Modified Julian Date for 0 hrs
]
    _sofa.iauCal2jd.restypes = c_int
except AttributeError:
    pass
cal2jd_msg = {0: 'OK', # Unused
                -1:'bad year',
                -2:'bad month',
                -3:'bad day'}
def cal2jd(iy, im, idd):
    """
    **  - - - - - - - - - -
    **   i a u C a l 2 j d
    **  - - - - - - - - - -
    **
    **  Gregorian Calendar to Julian Date.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 1)
    **
    **  Returned:
    **     djm0      double  MJD zero-point: always 2400000.5
    **     djm       double  Modified Julian Date for 0 hrs
    **
    **  Returned (function value):
    **               int     status:
    **                           0 = OK
    **                          -1 = bad year   (Note 3: JD not computed)
    **                          -2 = bad month  (JD not computed)
    **                          -3 = bad day    (JD computed)
    **
    **  Notes:
    **
    **  1) The algorithm used is valid from -4800 March 1, but this
    **     implementation rejects dates before -4799 January 1.
    **
    **  2) The Julian Date is returned in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding djm0 and
    **     djm.
    **
    **  3) In early eras the conversion is from the "Proleptic Gregorian
    **     Calendar";  no account is taken of the date(s) of adoption of
    **     the Gregorian Calendar, nor is the AD/BC numbering convention
    **     observed.
    **
    **  Reference:
    **
    **     Explanatory Supplement to the Astronomical Almanac,
    **     P. Kenneth Seidelmann (ed), University Science Books (1992),
    **     Section 12.92 (p604).
    **
    **  This revision:  2013 August 7
    """
    djm0 = c_double()
    djm = c_double()
    s=_sofa.iauCal2jd(iy, im, idd, byref(djm0), byref(djm))
    if s != 0:
        warnings.warn(cal2jd_msg[s], UserWarning, 2)
    return djm0.value, djm.value


#iauEpb
def epb (dj1,dj2):
    """
    *
    ** − − − − − − −
    ** i a u E p b
    ** − − − − − − −
    **
    ** Julian Date to Besselian Epoch.
    **
    ** This function is part of the International Astronomical Union’s
    ** SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    ** Status: support function.
    **
    ** Given:
    ** dj1,dj2 double Julian Date (see note)
    **
    ** Returned (function value):
    ** double Besselian Epoch.
    **
    ** Note:
    **
    ** The Julian Date is supplied in two pieces, in the usual SOFA
    ** manner, which is designed to preserve time resolution. The
    ** Julian Date is available as a single number by adding dj1 and
    ** dj2. The maximum resolution is achieved if dj1 is 2451545.0
    ** (J2000.0).
    **
    ** Reference:
    **
    ** Lieske, J.H., 1979. Astron.Astrophys., 73, 282

    """
    return 1900.0 + ((dj1 - DJ00) + (dj2 + 36524.68648)) / DTY
    


#iauEpb2jd
try:
    _sofa.iauEpb2jd.argtypes = [
    c_double,  # epb \Besselian Epoch (e.g. 1957.3)
    POINTER(c_double),  # djm0\ MJD zero−point: always 2400000.5
    POINTER(c_double),  # djm \ Modified Julian Date for 0 hrs
]
except AttributeError:
    pass
def epb2jd (epb):
    """
    ** − − − − − − − − − −
    ** i a u E p b 2 j d
    ** − − − − − − − − − −
    **
    ** Besselian Epoch to Julian Date.
    **
    ** This function is part of the International Astronomical Union’s
    ** SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    ** Status: support function.
    **
    ** Given:
    ** epb double Besselian Epoch (e.g. 1957.3)
    **
    ** Returned:
    ** djm0 double MJD zero−point: always 2400000.5
    ** djm double Modified Julian Date
    **
    ** Note:
    **
    ** The Julian Date is returned in two pieces, in the usual SOFA
    ** manner, which is designed to preserve time resolution. The
    ** Julian Date is available as a single number by adding djm0 and
    ** djm.
    **
    ** Reference:
    **
    ** Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
    """
    djm0 = c_double()
    djm = c_double()
    _sofa.iauEpb2jd(epb,byref(djm0), byref(djm))
    return djm0.value, djm.value


#iauEpj
def epj (dj1,dj2):
    """
    **  - - - - - - -
    **   i a u E p j
    **  - - - - - - -
    **
    **  Julian Date to Julian Epoch.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     dj1,dj2    double     Julian Date (see note)
    **
    **  Returned (function value):
    **                double     Julian Epoch
    **
    **  Note:
    **
    **     The Julian Date is supplied in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding dj1 and
    **     dj2.  The maximum resolution is achieved if dj1 is 2451545.0
    **     (J2000.0).
    **
    **  Reference:
    **
    **     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
    """
    return 2000.0 + ((dj1 - DJ00) + dj2) / DJY


#iauEpj2jd
try:
    _sofa.iauEpj2jd.argtypes = [
    c_double,  # epj \Julian Epoch (e.g. 1957.3)
    POINTER(c_double),  # djm0\ MJD zero−point: always 2400000.5
    POINTER(c_double),  # djm \ Modified Julian Date for 0 hrs
]
except AttributeError:
    pass
def epj2jd (epb):
    """
    **  - - - - - - - - - -
    **   i a u E p j 2 j d
    **  - - - - - - - - - -
    **
    **  Julian Epoch to Julian Date.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     epj      double    Julian Epoch (e.g. 1996.8)
    **
    **  Returned:
    **     djm0     double    MJD zero-point: always 2400000.5
    **     djm      double    Modified Julian Date
    **
    **  Note:
    **
    **     The Julian Date is returned in two pieces, in the usual SOFA
    **     manner, which is designed to preserve time resolution.  The
    **     Julian Date is available as a single number by adding djm0 and
    **     djm.
    **
    **  Reference:
    **
    **     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
    **
    **  This revision:  2013 August 7
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    djm0 = c_double()
    djm = c_double()
    _sofa.iauEpj2jd(epb,byref(djm0), byref(djm))
    return djm0.value, djm.value

#iauJd2cal
try:
    _sofa.iauJd2cal.argtypes = [
    c_double,  # jd1 \Julian Date (Notes 1, 2)
    c_double,  # jd2 \Julian Date (Notes 1, 2)
    POINTER(c_int), # iy \year
    POINTER(c_int),  # im \ month
    POINTER(c_int),  # idd \ day
    POINTER(c_double),  # fd \ fraction of day
]
    _sofa.iauJd2cal.restypes = c_int
except AttributeError:
    pass
jd2cal_msg = {0: 'OK', # Unused
             -1:' unacceptable date (Note 1)'}
def jd2cal(jd1, jd2):
    """
    ** − − − − − − − − − −
    ** i a u J d 2 c a l
    ** − − − − − − − − − −
    **
    ** Julian Date to Gregorian year, month, day, and fraction of a day.
    **
    ** This function is part of the International Astronomical Union’s
    ** SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    ** Status: support function.
    **
    ** Given:
    ** dj1,dj2 double Julian Date (Notes 1, 2)
    **
    ** Returned (arguments):
    ** iy int year
    ** im int month
    ** id int day
    ** fd double fraction of day
    **
    ** Returned (function value):
    ** int status:
    ** 0 = OK
    ** −1 = unacceptable date (Note 1)
    **
    ** Notes:
    **
    ** 1) The earliest valid date is −68569.5 (−4900 March 1). The
    ** largest value accepted is 1e9.
    **
    ** 2) The Julian Date is apportioned in any convenient way between
    ** the arguments dj1 and dj2. For example, JD=2450123.7 could
    ** be expressed in any of these ways, among others:
    **
    ** dj1 dj2
    **
    ** 2450123.7 0.0 (JD method)
    ** 2451545.0 −1421.3 (J2000 method)
    ** 2400000.5 50123.2 (MJD method)
    ** 2450123.5 0.2 (date & time method)
    **
    ** 3) In early eras the conversion is from the "proleptic Gregorian
    ** calendar"; no account is taken of the date(s) of adoption of
    ** the Gregorian calendar, nor is the AD/BC numbering convention
    ** observed.
    **
    ** Reference:
    **
    ** Explanatory Supplement to the Astronomical Almanac,
    ** P. Kenneth Seidelmann (ed), University Science Books (1992),
    """
    iy = c_int()
    im = c_int()
    idd = c_int()
    fd = c_double()
    
    s=_sofa.iauJd2cal(jd1, jd2, byref(iy), byref(im), byref(idd), byref(fd))
    if s != 0:
        warnings.warn(jd2cal_msg[s], UserWarning, 2)
    return iy.value, im.value, idd.value, fd.value

# iauJdcalf
jdcalf_msg = {0: 'OK', # Unused
              -1:' date out of range',
              +1:'NDP not 0−9 (interpreted as 0)'}
try:
    _sofa.iauJdcalf.argtypes = [
    c_int, # ndp \number of decimal places of days in fraction
    c_double,  # jd1 \Julian Date (Notes 1, 2)
    c_double,  # jd2 \Julian Date (Notes 1, 2)
    c_int * 4 #iymdf \year, month, day, fraction in Gregorian calendar
    
]
    _sofa.iauJdcalf.restypes = c_int
except AttributeError:
    pass  
def jdcalf(ndp,jd1,jd2):
    """
    

    **  - - - - - - - - - -
    **   i a u J d c a l f
    **  - - - - - - - - - -
    **
    **  Julian Date to Gregorian Calendar, expressed in a form convenient
    **  for formatting messages:  rounded to a specified precision.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ndp       int      number of decimal places of days in fraction
    **     dj1,dj2   double   dj1+dj2 = Julian Date (Note 1)
    **
    **  Returned:
    **     iymdf     int[4]   year, month, day, fraction in Gregorian
    **                        calendar
    **
    **  Returned (function value):
    **               int      status:
    **                          -1 = date out of range
    **                           0 = OK
    **                          +1 = NDP not 0-9 (interpreted as 0)
    **
    **  Notes:
    **
    **  1) The Julian Date is apportioned in any convenient way between
    **     the arguments dj1 and dj2.  For example, JD=2450123.7 could
    **     be expressed in any of these ways, among others:
    **
    **             dj1            dj2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **  2) In early eras the conversion is from the "Proleptic Gregorian
    **     Calendar";  no account is taken of the date(s) of adoption of
    **     the Gregorian Calendar, nor is the AD/BC numbering convention
    **     observed.
    **
    **  3) Refer to the function iauJd2cal.
    **
    **  4) NDP should be 4 or less if internal overflows are to be
    **     avoided on machines which use 16-bit integers.
    **
    **  Called:
    **     iauJd2cal    JD to Gregorian calendar
    **
    **  Reference:
    **
    **     Explanatory Supplement to the Astronomical Almanac,
    **     P. Kenneth Seidelmann (ed), University Science Books (1992),
    **     Section 12.92 (p604).

    """
    iymdf = (c_int * 4)()
    s = _sofa.iauJdcalf(ndp, jd1, jd2, iymdf)
    if s != 0:
        warnings.warn(jdcalf_msg[s], UserWarning, 2)
    return tuple([v for v in iymdf])



###########################################################################       Astrometry
################################################################

# iauAb
try:
    _sofa.iauAb.argtypes = [
                            ndpointer(shape=3, dtype=float), #pnat /natural direction to the source (unit vector)
                            ndpointer(shape=3, dtype=float), #v /observer barycentric velocity in units of c
                            c_double, #s /distance between the Sun and the observer (au)
                            c_double, #bml /sqrt(1−|v|^2): reciprocal of Lorenz factor
                            ndpointer(shape=3, dtype=float)]#ppr /proper direction to source (unit vector)
except AttributeError:
    pass  
def ab(pnat,v,s,bml):
    """

    ** − − − − − −
    ** i a u A b
    ** − − − − − −
    **
    ** Apply aberration to transform natural direction into proper
    ** direction.
    **
    ** This function is part of the International Astronomical Union’s
    ** SOFA (Standards of Fundamental Astronomy) software collection.
    **
    ** Status: support function.
    **
    ** Given:
    ** pnat double[3] natural direction to the source (unit vector)
    ** v double[3] observer barycentric velocity in units of c
    ** s double distance between the Sun and the observer (au)
    ** bm1 double sqrt(1−|v|^2): reciprocal of Lorenz factor
    **
    ** Returned:
    ** ppr double[3] proper direction to source (unit vector)
    **
    ** Notes:
    **
    ** 1) The algorithm is based on Expr. (7.40) in the Explanatory
    ** Supplement (Urban & Seidelmann 2013), but with the following
    ** changes:
    **
    ** o Rigorous rather than approximate normalization is applied.
    **
    ** o The gravitational potential term from Expr. (7) in
    ** Klioner (2003) is added, taking into account only the Sun’s
    ** contribution. This has a maximum effect of about
    ** 0.4 microarcsecond.
    **
    ** 2) In almost all cases, the maximum accuracy will be limited by the
    ** supplied velocity. For example, if the SOFA iauEpv00 function is
    ** used, errors of up to 5 microarcseconds could occur.
    **
    ** References:
    **
    ** Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    ** the Astronomical Almanac, 3rd ed., University Science Books
    ** (2013).
    **
    ** Klioner, Sergei A., "A practical relativistic model for micro−
    ** arcsecond astrometry in space", Astr. J. 125, 1580−1597 (2003).

    """
    ppr = zeros(shape=3, dtype=float)
    _sofa.iauAb(pnat,v,s,bml,ppr)
    return ppr

# iauApcg
try:
    _sofa.iauApcg.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1)
                            ndpointer(shape=(2,3), dtype=float), #ebpv /Earth barycentric pos/vel (au, au/day)
                            ndpointer(shape=3, dtype=float), #ehp /Earth heliocentric position (au)
                            POINTER(ASTROM)]#astrom Star-independent astrometry parameters 
except AttributeError:
    pass  
def apcg(date1,date2,ebpv,ehp):
    """
    **  - - - - - - - -
    **   i a u A p c g
    **  - - - - - - - -
    **
    **  For a geocentric observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and GCRS coordinates.
    **  The Earth ephemeris is supplied by the caller.
    **
    **  The parameters produced by this function are required in the
    **  parallax, light deflection and aberration parts of the astrometric
    **  transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double       TDB as a 2-part...
    **     date2  double       ...Julian Date (Note 1)
    **     ebpv   double[2][3] Earth barycentric pos/vel (au, au/day)
    **     ehp    double[3]    Earth heliocentric position (au)
    **
    **  Returned:
    **     astrom iauASTROM*   star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  4) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauApcs      astrometry parameters, ICRS-GCRS, space observer
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22


    """
    astrom  =ASTROM()
    _sofa.iauApcg(date1,date2,ebpv,ehp,astrom)
    return astrom

#iauApcg13
try:
    _sofa.iauApcg13.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1)
                            POINTER(ASTROM)]#astrom Star-independent astrometry parameters 
except AttributeError:
    pass  
def apcg13(date1,date2):
    """
    **  - - - - - - - - - -
    **   i a u A p c g 1 3
    **  - - - - - - - - - -
    **
    **  For a geocentric observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and GCRS coordinates.
    **  The caller supplies the date, and SOFA models are used to predict
    **  the Earth ephemeris.
    **
    **  The parameters produced by this function are required in the
    **  parallax, light deflection and aberration parts of the astrometric
    **  transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double     TDB as a 2-part...
    **     date2  double     ...Julian Date (Note 1)
    **
    **  Returned:
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) In cases where the caller wishes to supply his own Earth
    **     ephemeris, the function iauApcg can be used instead of the present
    **     function.
    **
    **  4) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  5) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauEpv00     Earth position and velocity
    **     iauApcg      astrometry parameters, ICRS-GCRS, geocenter

    """
    astrom  =ASTROM()
    _sofa.iauApcg13(date1,date2,astrom)
    return astrom


#iauApci
try:
    _sofa.iauApci.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1) 
                            ndpointer(shape=(2,3), dtype=float), #ebpv /Earth barycentric position/velocity (au, au/day)
                            ndpointer(shape=3, dtype=float),#ehp /Earth heliocentric position (au)
                            c_double,#x, /CIP X,Y (components of unit vector)
                            c_double,#y, /CIP X,Y (components of unit vector)
                            c_double,#s /the CIO locator s (radians)
                            POINTER(ASTROM)]
except AttributeError:
    pass  
def apci(date1,date2,ebpv,ehp,x,y,s):
    """
    **  - - - - - - - -
    **   i a u A p c i
    **  - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and geocentric CIRS
    **  coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
    **  caller.
    **
    **  The parameters produced by this function are required in the
    **  parallax, light deflection, aberration, and bias-precession-nutation
    **  parts of the astrometric transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double       TDB as a 2-part...
    **     date2  double       ...Julian Date (Note 1)
    **     ebpv   double[2][3] Earth barycentric position/velocity (au, au/day)
    **     ehp    double[3]    Earth heliocentric position (au)
    **     x,y    double       CIP X,Y (components of unit vector)
    **     s      double       the CIO locator s (radians)
    **
    **  Returned:
    **     astrom iauASTROM*   star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) In cases where the caller does not wish to provide the Earth
    **     ephemeris and CIP/CIO, the function iauApci13 can be used instead
    **     of the present function.  This computes the required quantities
    **     using other SOFA functions.
    **
    **  4) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  5) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauApcg      astrometry parameters, ICRS-GCRS, geocenter
    **     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
    **
    **  This revision:   2013 September 25
    """
    astrom  =ASTROM()
    _sofa.iauApci(date1,date2,ebpv,ehp,x,y,s,astrom)
    return astrom


#iauApci13
try:
    _sofa.iauApci13.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1) 
                            POINTER(ASTROM),# astrom /star−independent astrometry parameters:
                            POINTER(c_double) #eo /equation of the origins (ERA−GST)
                            ]
except AttributeError:
    pass  
def apci13(date1,date2):
    """
    **  - - - - - - - - - -
    **   i a u A p c i 1 3
    **  - - - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and geocentric CIRS
    **  coordinates.  The caller supplies the date, and SOFA models are used
    **  to predict the Earth ephemeris and CIP/CIO.
    **
    **  The parameters produced by this function are required in the
    **  parallax, light deflection, aberration, and bias-precession-nutation
    **  parts of the astrometric transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double      TDB as a 2-part...
    **     date2  double      ...Julian Date (Note 1)
    **
    **  Returned:
    **     astrom iauASTROM*  star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **     eo     double*     equation of the origins (ERA-GST)
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) In cases where the caller wishes to supply his own Earth
    **     ephemeris and CIP/CIO, the function iauApci can be used instead
    **     of the present function.
    **
    **  4) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  5) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauEpv00     Earth position and velocity
    **     iauPnm06a    classical NPB matrix, IAU 2006/2000A
    **     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
    **     iauS06       the CIO locator s, given X,Y, IAU 2006
    **     iauApci      astrometry parameters, ICRS-CIRS
    **     iauEors      equation of the origins, given NPB matrix and s
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    astrom  =ASTROM()
    eo = c_double()
    _sofa.iauApci13(date1, date2, astrom, byref(eo))
    return astrom, eo.value

# iauApco
try:
    _sofa.iauApco.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1) 
                            ndpointer(shape=(2,3), dtype=float), #ebpv /Earth barycentric position/velocity (au, au/day)
                            ndpointer(shape=3, dtype=float),#ehp /Earth heliocentric position (au)
                            c_double,#x, /CIP X,Y (components of unit vector)
                            c_double,#y, /CIP X,Y (components of unit vector)
                            c_double,#s /the CIO locator s (radians)
                            c_double,# theta /Earth rotation angle (radians)
                            c_double, # elong /longitude (radians, east +ve, Note 3)
                            c_double, # phi /latitude (geodetic, radians, Note 3)
                            c_double, # hm / height above ellipsoid (m, geodetic, Note 3)
                            c_double, # xp /polar motion coordinates (radians, Note 4)
                            c_double, # yp /polar motion coordinates (radians, Note 4)
                            c_double, # sp /the TIO locator s' (radians, Note 4)
                            c_double, # refa /refraction constant A (radians, Note 5)
                            c_double, # refb /refraction constant B (radians, Note 5)
                            POINTER(ASTROM)]
except AttributeError:
    pass  
def apco(date1,date2,ebpv,ehp,x,y,s,theta,elong,phi,hm,xp,yp, sp, refa, refb):
    """
    **  - - - - - - - -
    **   i a u A p c o
    **  - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and observed
    **  coordinates.  The caller supplies the Earth ephemeris, the Earth
    **  rotation information and the refraction constants as well as the
    **  site coordinates.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double       TDB as a 2-part...
    **     date2  double       ...Julian Date (Note 1)
    **     ebpv   double[2][3] Earth barycentric PV (au, au/day, Note 2)
    **     ehp    double[3]    Earth heliocentric P (au, Note 2)
    **     x,y    double       CIP X,Y (components of unit vector)
    **     s      double       the CIO locator s (radians)
    **     theta  double       Earth rotation angle (radians)
    **     elong  double       longitude (radians, east +ve, Note 3)
    **     phi    double       latitude (geodetic, radians, Note 3)
    **     hm     double       height above ellipsoid (m, geodetic, Note 3)
    **     xp,yp  double       polar motion coordinates (radians, Note 4)
    **     sp     double       the TIO locator s' (radians, Note 4)
    **     refa   double       refraction constant A (radians, Note 5)
    **     refb   double       refraction constant B (radians, Note 5)
    **
    **  Returned:
    **     astrom iauASTROM*   star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) The vectors eb, eh, and all the astrom vectors, are with respect
    **     to BCRS axes.
    **
    **  3) The geographical coordinates are with respect to the WGS84
    **     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
    **     CONVENTION:  the longitude required by the present function is
    **     right-handed, i.e. east-positive, in accordance with geographical
    **     convention.
    **
    **  4) xp and yp are the coordinates (in radians) of the Celestial
    **     Intermediate Pole with respect to the International Terrestrial
    **     Reference System (see IERS Conventions), measured along the
    **     meridians 0 and 90 deg west respectively.  sp is the TIO locator
    **     s', in radians, which positions the Terrestrial Intermediate
    **     Origin on the equator.  For many applications, xp, yp and
    **     (especially) sp can be set to zero.
    **
    **     Internally, the polar motion is stored in a form rotated onto the
    **     local meridian.
    **
    **  5) The refraction constants refa and refb are for use in a
    **     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
    **     (i.e. refracted) zenith distance and dZ is the amount of
    **     refraction.
    **
    **  6) It is advisable to take great care with units, as even unlikely
    **     values of the input parameters are accepted and processed in
    **     accordance with the models used.
    **
    **  7) In cases where the caller does not wish to provide the Earth
    **     Ephemeris, the Earth rotation information and refraction
    **     constants, the function iauApco13 can be used instead of the
    **     present function.  This starts from UTC and weather readings etc.
    **     and computes suitable values using other SOFA functions.
    **
    **  8) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  9) The context structure astrom produced by this function is used by
    **     iauAtioq, iauAtoiq, iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauAper      astrometry parameters: update ERA
    **     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
    **     iauPvtob     position/velocity of terrestrial station
    **     iauTrxpv     product of transpose of r-matrix and pv-vector
    **     iauApcs      astrometry parameters, ICRS-GCRS, space observer
    **     iauCr        copy r-matrix
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    astrom  =ASTROM()
    _sofa.iauApco(date1,date2,ebpv,ehp,x,y,s,theta,elong,phi,hm,xp,yp, sp, refa, refb,astrom)
    return astrom


# iauApco13
try:
    _sofa.iauApco13.argtypes = [
                            c_double, #utc1 /TDB as a 2−part
                            c_double, #utc2 /...quasi Julian Date (Notes 1,2)
                            c_double, #dut1 /UT1-UTC (seconds, Note 3)
                            c_double, # elong /longitude (radians, east +ve, Note 3)
                            c_double, # phi /latitude (geodetic, radians, Note 3)
                            c_double, # hm / height above ellipsoid (m, geodetic, Note 3)
                            c_double, # xp /polar motion coordinates (radians, Note 4)
                            c_double, # yp /polar motion coordinates (radians, Note 4)
                            c_double, # phpa /pressure at the observer (hPa = mB, Note 6)
                            c_double, # tc /ambient temperature at the observer (deg C)
                            c_double, # rh / relative humidity at the observer (range 0-1)
                            c_double, # wl / wavelength (micrometers, Note 7)
                            POINTER(ASTROM),# astrom /star-independent astrometry parameters:
                            POINTER(c_double)]# eo /equation of the origins (ERA-GST)
except AttributeError:
    pass  
apco13_msg = {0: 'OK', # Unused
                +1:'dubious year (Note 2)',
                -1:'unacceptable date'}
def apco13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tc,rh,wl):
    """
    **  - - - - - - - - - -
    **   i a u A p c o 1 3
    **  - - - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between ICRS and observed
    **  coordinates.  The caller supplies UTC, site coordinates, ambient air
    **  conditions and observing wavelength, and SOFA models are used to
    **  obtain the Earth ephemeris, CIP/CIO and refraction constants.
    **
    **  The parameters produced by this function are required in the
    **  parallax, light deflection, aberration, and bias-precession-nutation
    **  parts of the ICRS/CIRS transformations.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     utc1   double     UTC as a 2-part...
    **     utc2   double     ...quasi Julian Date (Notes 1,2)
    **     dut1   double     UT1-UTC (seconds, Note 3)
    **     elong  double     longitude (radians, east +ve, Note 4)
    **     phi    double     latitude (geodetic, radians, Note 4)
    **     hm     double     height above ellipsoid (m, geodetic, Notes 4,6)
    **     xp,yp  double     polar motion coordinates (radians, Note 5)
    **     phpa   double     pressure at the observer (hPa = mB, Note 6)
    **     tc     double     ambient temperature at the observer (deg C)
    **     rh     double     relative humidity at the observer (range 0-1)
    **     wl     double     wavelength (micrometers, Note 7)
    **
    **  Returned:
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **     eo     double*    equation of the origins (ERA-GST)
    **
    **  Returned (function value):
    **            int        status: +1 = dubious year (Note 2)
    **                                0 = OK
    **                               -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  2)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the
    **      future to be trusted.  See iauDat for further details.
    **
    **  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  4)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many
    **      applications, xp and yp can be set to zero.
    **
    **      Internally, the polar motion is stored in a form rotated onto
    **      the local meridian.
    **
    **  6)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB), is
    **      available, an adequate estimate of hm can be obtained from the
    **      expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to
    **      the pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  7)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  8)  It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  9)  In cases where the caller wishes to supply his own Earth
    **      ephemeris, Earth rotation information and refraction constants,
    **      the function iauApco can be used instead of the present function.
    **
    **  10) This is one of several functions that inserts into the astrom
    **      structure star-independent parameters needed for the chain of
    **      astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **      The various functions support different classes of observer and
    **      portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **      Those with names ending in "13" use contemporary SOFA models to
    **      compute the various ephemerides.  The others accept ephemerides
    **      supplied by the caller.
    **
    **      The transformation from ICRS to GCRS covers space motion,
    **      parallax, light deflection, and aberration.  From GCRS to CIRS
    **      comprises frame bias and precession-nutation.  From CIRS to
    **      observed takes account of Earth rotation, polar motion, diurnal
    **      aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **      transformation), and atmospheric refraction.
    **
    **  11) The context structure astrom produced by this function is used
    **      by iauAtioq, iauAtoiq, iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauUtctai    UTC to TAI
    **     iauTaitt     TAI to TT
    **     iauUtcut1    UTC to UT1
    **     iauEpv00     Earth position and velocity
    **     iauPnm06a    classical NPB matrix, IAU 2006/2000A
    **     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
    **     iauS06       the CIO locator s, given X,Y, IAU 2006
    **     iauEra00     Earth rotation angle, IAU 2000
    **     iauSp00      the TIO locator s', IERS 2000
    **     iauRefco     refraction constants for given ambient conditions
    **     iauApco      astrometry parameters, ICRS-observed
    **     iauEors      equation of the origins, given NPB matrix and s
    **
    **  This revision:   2013 December 5
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    astrom  =ASTROM()
    eo = c_double()
    s = _sofa.iauApco13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tc,rh,wl,astrom,byref(eo))
    if s != 0:
        warnings.warn(apco13_msg[s], UserWarning, 2)
    return astrom, eo.value

#iauApcs
try:
    _sofa.iauApcs.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1) 
                            ndpointer(shape=(2,3), dtype=float),# pv /observer's geocentric pos/vel (m, m/s)
                            ndpointer(shape=(2,3), dtype=float), #ebpv /Earth barycentric position/velocity (au, au/day)
                            ndpointer(shape=3, dtype=float),#ehp /Earth heliocentric position (au)
                            POINTER(ASTROM)]
except AttributeError:
    pass  
def apcs(date1,date2,pv,ebpv,ehp):
    """
    **  - - - - - - - -
    **   i a u A p c s
    **  - - - - - - - -
    **
    **  For an observer whose geocentric position and velocity are known,
    **  prepare star-independent astrometry parameters for transformations
    **  between ICRS and GCRS.  The Earth ephemeris is supplied by the
    **  caller.
    **
    **  The parameters produced by this function are required in the space
    **  motion, parallax, light deflection and aberration parts of the
    **  astrometric transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double       TDB as a 2-part...
    **     date2  double       ...Julian Date (Note 1)
    **     pv     double[2][3] observer's geocentric pos/vel (m, m/s)
    **     ebpv   double[2][3] Earth barycentric PV (au, au/day)
    **     ehp    double[3]    Earth heliocentric P (au)
    **
    **  Returned:
    **     astrom iauASTROM*   star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) Providing separate arguments for (i) the observer's geocentric
    **     position and velocity and (ii) the Earth ephemeris is done for
    **     convenience in the geocentric, terrestrial and Earth orbit cases.
    **     For deep space applications it maybe more convenient to specify
    **     zero geocentric position and velocity and to supply the
    **     observer's position and velocity information directly instead of
    **     with respect to the Earth.  However, note the different units:
    **     m and m/s for the geocentric vectors, au and au/day for the
    **     heliocentric and barycentric vectors.
    **
    **  4) In cases where the caller does not wish to provide the Earth
    **     ephemeris, the function iauApcs13 can be used instead of the
    **     present function.  This computes the Earth ephemeris using the
    **     SOFA function iauEpv00.
    **
    **  5) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  6) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauCp        copy p-vector
    **     iauPm        modulus of p-vector
    **     iauPn        decompose p-vector into modulus and direction
    **     iauIr        initialize r-matrix to identity
    **
    **  This revision:   2017 March 16
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    astrom = ASTROM()
    _sofa.iauApcs(date1,date2,pv,ebpv,ehp,astrom)
    return astrom

#iauApcs13
try:
    _sofa.iauApcs13.argtypes = [
                            c_double, #date1 /TDB as a 2−part
                            c_double, #date2 /Julian Date (Note 1) 
                            ndpointer(shape=(2,3), dtype=float),# pv /observer's geocentric pos/vel (m, m/s)
                            POINTER(ASTROM)] # astrom, /star-independent astrometry parameters:
except AttributeError:
    pass  
def apcs13(date1,date2,pv):
    """
    **  - - - - - - - - - -
    **   i a u A p c s 1 3
    **  - - - - - - - - - -
    **
    **  For an observer whose geocentric position and velocity are known,
    **  prepare star-independent astrometry parameters for transformations
    **  between ICRS and GCRS.  The Earth ephemeris is from SOFA models.
    **
    **  The parameters produced by this function are required in the space
    **  motion, parallax, light deflection and aberration parts of the
    **  astrometric transformation chain.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1  double       TDB as a 2-part...
    **     date2  double       ...Julian Date (Note 1)
    **     pv     double[2][3] observer's geocentric pos/vel (Note 3)
    **
    **  Returned:
    **     astrom iauASTROM*   star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       unchanged
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) All the vectors are with respect to BCRS axes.
    **
    **  3) The observer's position and velocity pv are geocentric but with
    **     respect to BCRS axes, and in units of m and m/s.  No assumptions
    **     are made about proximity to the Earth, and the function can be
    **     used for deep space applications as well as Earth orbit and
    **     terrestrial.
    **
    **  4) In cases where the caller wishes to supply his own Earth
    **     ephemeris, the function iauApcs can be used instead of the present
    **     function.
    **
    **  5) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  6) The context structure astrom produced by this function is used by
    **     iauAtciq* and iauAticq*.
    **
    **  Called:
    **     iauEpv00     Earth position and velocity
    **     iauApcs      astrometry parameters, ICRS-GCRS, space observer
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    {
       double ehpv[2][3], ebpv[2][3];
    
    
    /* Earth barycentric & heliocentric position/velocity (au, au/d). */
       (void) iauEpv00(date1, date2, ehpv, ebpv);
    
    /* Compute the star-independent astrometry parameters. */
       iauApcs(date1, date2, pv, ebpv, ehpv[0], astrom);
    
    /* Finished. */
    """
    astrom = ASTROM()
    _sofa.iauApcs13(date1,date2,pv,astrom)
    return astrom

#iauAper
try:
    _sofa.iauAper.argtypes = [
                            c_double, #theta /Earth rotation angle (radians, Note 2)
                            POINTER(ASTROM)] # astrom, /star-independent astrometry parameters:
except AttributeError:
    pass  
def aper(theta,astrom):
    """
    /*
    **  - - - - - - - -
    **   i a u A p e r
    **  - - - - - - - -
    **
    **  In the star-independent astrometry parameters, update only the
    **  Earth rotation angle, supplied by the caller explicitly.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     theta   double      Earth rotation angle (radians, Note 2)
    **     astrom  iauASTROM*  star-independent astrometry parameters:
    **      pmt    double       not used
    **      eb     double[3]    not used
    **      eh     double[3]    not used
    **      em     double       not used
    **      v      double[3]    not used
    **      bm1    double       not used
    **      bpn    double[3][3] not used
    **      along  double       longitude + s' (radians)
    **      xpl    double       not used
    **      ypl    double       not used
    **      sphi   double       not used
    **      cphi   double       not used
    **      diurab double       not used
    **      eral   double       not used
    **      refa   double       not used
    **      refb   double       not used
    **
    **  Returned:
    **     astrom  iauASTROM*  star-independent astrometry parameters:
    **      pmt    double       unchanged
    **      eb     double[3]    unchanged
    **      eh     double[3]    unchanged
    **      em     double       unchanged
    **      v      double[3]    unchanged
    **      bm1    double       unchanged
    **      bpn    double[3][3] unchanged
    **      along  double       unchanged
    **      xpl    double       unchanged
    **      ypl    double       unchanged
    **      sphi   double       unchanged
    **      cphi   double       unchanged
    **      diurab double       unchanged
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       unchanged
    **      refb   double       unchanged
    **
    **  Notes:
    **
    **  1) This function exists to enable sidereal-tracking applications to
    **     avoid wasteful recomputation of the bulk of the astrometry
    **     parameters:  only the Earth rotation is updated.
    **
    **  2) For targets expressed as equinox based positions, such as
    **     classical geocentric apparent (RA,Dec), the supplied theta can be
    **     Greenwich apparent sidereal time rather than Earth rotation
    **     angle.
    **
    **  3) The function iauAper13 can be used instead of the present
    **     function, and starts from UT1 rather than ERA itself.
    **
    **  4) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  This revision:   2013 September 25
    **
    **  SOFA release 2019-07-22
    """
    _sofa.iauAper(theta,astrom)
    return astrom

# iauAper13

_sofa.iauAper13.argtypes = [
                            c_double, #utt1 /TDB as a 2−part
                            c_double, #ut12 /Julian Date (Note 1) 
                            POINTER(ASTROM)] # astrom, /star-independent astrometry parameters:

def aper13(ut11,ut12,astrom):
    """
    **  - - - - - - - - - -
**   i a u A p c s 1 3
**  - - - - - - - - - -
**
**  For an observer whose geocentric position and velocity are known,
**  prepare star-independent astrometry parameters for transformations
**  between ICRS and GCRS.  The Earth ephemeris is from SOFA models.
**
**  The parameters produced by this function are required in the space
**  motion, parallax, light deflection and aberration parts of the
**  astrometric transformation chain.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1  double       TDB as a 2-part...
**     date2  double       ...Julian Date (Note 1)
**     pv     double[2][3] observer's geocentric pos/vel (Note 3)
**
**  Returned:
**     astrom iauASTROM*   star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       unchanged
**      xpl    double       unchanged
**      ypl    double       unchanged
**      sphi   double       unchanged
**      cphi   double       unchanged
**      diurab double       unchanged
**      eral   double       unchanged
**      refa   double       unchanged
**      refb   double       unchanged
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways, among
**     others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  For most
**     applications of this function the choice will not be at all
**     critical.
**
**     TT can be used instead of TDB without any significant impact on
**     accuracy.
**
**  2) All the vectors are with respect to BCRS axes.
**
**  3) The observer's position and velocity pv are geocentric but with
**     respect to BCRS axes, and in units of m and m/s.  No assumptions
**     are made about proximity to the Earth, and the function can be
**     used for deep space applications as well as Earth orbit and
**     terrestrial.
**
**  4) In cases where the caller wishes to supply his own Earth
**     ephemeris, the function iauApcs can be used instead of the present
**     function.
**
**  5) This is one of several functions that inserts into the astrom
**     structure star-independent parameters needed for the chain of
**     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
**
**     The various functions support different classes of observer and
**     portions of the transformation chain:
**
**          functions         observer        transformation
**
**       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
**       iauApci iauApci13    terrestrial     ICRS <-> CIRS
**       iauApco iauApco13    terrestrial     ICRS <-> observed
**       iauApcs iauApcs13    space           ICRS <-> GCRS
**       iauAper iauAper13    terrestrial     update Earth rotation
**       iauApio iauApio13    terrestrial     CIRS <-> observed
**
**     Those with names ending in "13" use contemporary SOFA models to
**     compute the various ephemerides.  The others accept ephemerides
**     supplied by the caller.
**
**     The transformation from ICRS to GCRS covers space motion,
**     parallax, light deflection, and aberration.  From GCRS to CIRS
**     comprises frame bias and precession-nutation.  From CIRS to
**     observed takes account of Earth rotation, polar motion, diurnal
**     aberration and parallax (unless subsumed into the ICRS <-> GCRS
**     transformation), and atmospheric refraction.
**
**  6) The context structure astrom produced by this function is used by
**     iauAtciq* and iauAticq*.
**
**  Called:
**     iauEpv00     Earth position and velocity
**     iauApcs      astrometry parameters, ICRS-GCRS, space observer
**
**  This revision:   2013 October 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    _sofa.iauAper13(ut11,ut12,astrom)
    return astrom

#iauApio
try:
    _sofa.iauApio.argtypes = [
                            c_double, # sp /the TIO locator s' (radians, Note 1)
                            c_double, # theta /Earth rotation angle (radians) 
                            c_double, # elong /longitude (radians, east +ve, Note 2)
                            c_double, # phi /geodetic latitude (radians, Note 2)
                            c_double, # hm /height above ellipsoid (m, geodetic Note 2)
                            c_double, # xp /polar motion coordinates (radians, Note 3)
                            c_double, # yp /polar motion coordinates (radians, Note 3)
                            c_double, # refa /refraction constant A (radians, Note 4)
                            c_double, # refb /refraction constant B (radians, Note 4)
                            POINTER(ASTROM)] # astrom, /star-independent astrometry parameters:
except AttributeError:
    pass  
def apio(sp,theta,elong,phi,hm,xp,yp,refa,refb):
    """
    **  - - - - - - - -
    **   i a u A p i o
    **  - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between CIRS and observed
    **  coordinates.  The caller supplies the Earth orientation information
    **  and the refraction constants as well as the site coordinates.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     sp     double      the TIO locator s' (radians, Note 1)
    **     theta  double      Earth rotation angle (radians)
    **     elong  double      longitude (radians, east +ve, Note 2)
    **     phi    double      geodetic latitude (radians, Note 2)
    **     hm     double      height above ellipsoid (m, geodetic Note 2)
    **     xp,yp  double      polar motion coordinates (radians, Note 3)
    **     refa   double      refraction constant A (radians, Note 4)
    **     refb   double      refraction constant B (radians, Note 4)
    **
    **  Returned:
    **     astrom iauASTROM*  star-independent astrometry parameters:
    **      pmt    double       unchanged
    **      eb     double[3]    unchanged
    **      eh     double[3]    unchanged
    **      em     double       unchanged
    **      v      double[3]    unchanged
    **      bm1    double       unchanged
    **      bpn    double[3][3] unchanged
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Notes:
    **
    **  1) sp, the TIO locator s', is a tiny quantity needed only by the
    **     most precise applications.  It can either be set to zero or
    **     predicted using the SOFA function iauSp00.
    **
    **  2) The geographical coordinates are with respect to the WGS84
    **     reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **     longitude required by the present function is east-positive
    **     (i.e. right-handed), in accordance with geographical convention.
    **
    **  3) The polar motion xp,yp can be obtained from IERS bulletins.  The
    **     values are the coordinates (in radians) of the Celestial
    **     Intermediate Pole with respect to the International Terrestrial
    **     Reference System (see IERS Conventions 2003), measured along the
    **     meridians 0 and 90 deg west respectively.  For many applications,
    **     xp and yp can be set to zero.
    **
    **     Internally, the polar motion is stored in a form rotated onto the
    **     local meridian.
    **
    **  4) The refraction constants refa and refb are for use in a
    **     dZ = A*tan(Z)+B*tan^3(Z) model, where Z is the observed
    **     (i.e. refracted) zenith distance and dZ is the amount of
    **     refraction.
    **
    **  5) It is advisable to take great care with units, as even unlikely
    **     values of the input parameters are accepted and processed in
    **     accordance with the models used.
    **
    **  6) In cases where the caller does not wish to provide the Earth
    **     rotation information and refraction constants, the function
    **     iauApio13 can be used instead of the present function.  This
    **     starts from UTC and weather readings etc. and computes suitable
    **     values using other SOFA functions.
    **
    **  7) This is one of several functions that inserts into the astrom
    **     structure star-independent parameters needed for the chain of
    **     astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **     The various functions support different classes of observer and
    **     portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **     Those with names ending in "13" use contemporary SOFA models to
    **     compute the various ephemerides.  The others accept ephemerides
    **     supplied by the caller.
    **
    **     The transformation from ICRS to GCRS covers space motion,
    **     parallax, light deflection, and aberration.  From GCRS to CIRS
    **     comprises frame bias and precession-nutation.  From CIRS to
    **     observed takes account of Earth rotation, polar motion, diurnal
    **     aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **     transformation), and atmospheric refraction.
    **
    **  8) The context structure astrom produced by this function is used by
    **     iauAtioq and iauAtoiq.
    **
    **  Called:
    **     iauPvtob     position/velocity of terrestrial station
    **     iauAper      astrometry parameters: update ERA
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    """
    astrom = ASTROM()
    _sofa.iauApio(sp,theta,elong,phi,hm,xp,yp,refa,refb,astrom)
    return astrom


#iauApio13
apio13_msg = {0: 'OK', # Unused
                +1:'dubious year (Note 2)',
                -1:'unacceptable date'}
try:
    _sofa.iauApio13.argtypes = [
                            c_double, # utc1 /UTC as a 2-part...
                            c_double, # utc2 /...quasi Julian Date (Notes 1,2) 
                            c_double, # dut1 /UT1-UTC (seconds)
                            c_double, # elong /longitude (radians, east +ve, Note 3)
                            c_double, # phi /geodetic latitude (radians, Note 3)
                            c_double, # hm /height above ellipsoid (m, geodetic Notes 4,6)
                            c_double, # xp /polar motion coordinates (radians, Note 5)
                            c_double, # yp /polar motion coordinates (radians, Note 3)
                            c_double, # phpa /pressure at the observer (hPa = mB, Note 6)
                            c_double, # tc /ambient temperature at the observer (deg C)
                            c_double, # rh /relative humidity at the observer (range 0-1)
                            c_double, # wl /wavelength (micrometers, Note 7)
                            POINTER(ASTROM)] # astrom, /star-independent astrometry parameters:
except AttributeError:
    pass  
def apio13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tc,rh,wl):
    """
    /*
    **  - - - - - - - - - -
    **   i a u A p i o 1 3
    **  - - - - - - - - - -
    **
    **  For a terrestrial observer, prepare star-independent astrometry
    **  parameters for transformations between CIRS and observed
    **  coordinates.  The caller supplies UTC, site coordinates, ambient air
    **  conditions and observing wavelength.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     utc1   double      UTC as a 2-part...
    **     utc2   double      ...quasi Julian Date (Notes 1,2)
    **     dut1   double      UT1-UTC (seconds)
    **     elong  double      longitude (radians, east +ve, Note 3)
    **     phi    double      geodetic latitude (radians, Note 3)
    **     hm     double      height above ellipsoid (m, geodetic Notes 4,6)
    **     xp,yp  double      polar motion coordinates (radians, Note 5)
    **     phpa   double      pressure at the observer (hPa = mB, Note 6)
    **     tc     double      ambient temperature at the observer (deg C)
    **     rh     double      relative humidity at the observer (range 0-1)
    **     wl     double      wavelength (micrometers, Note 7)
    **
    **  Returned:
    **     astrom iauASTROM*  star-independent astrometry parameters:
    **      pmt    double       unchanged
    **      eb     double[3]    unchanged
    **      eh     double[3]    unchanged
    **      em     double       unchanged
    **      v      double[3]    unchanged
    **      bm1    double       unchanged
    **      bpn    double[3][3] unchanged
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned (function value):
    **            int         status: +1 = dubious year (Note 2)
    **                                 0 = OK
    **                                -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  2)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the future
    **      to be trusted.  See iauDat for further details.
    **
    **  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  4)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many applications,
    **      xp and yp can be set to zero.
    **
    **      Internally, the polar motion is stored in a form rotated onto
    **      the local meridian.
    **
    **  6)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB), is
    **      available, an adequate estimate of hm can be obtained from the
    **      expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to the
    **      pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  7)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  8)  It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  9)  In cases where the caller wishes to supply his own Earth
    **      rotation information and refraction constants, the function
    **      iauApc can be used instead of the present function.
    **
    **  10) This is one of several functions that inserts into the astrom
    **      structure star-independent parameters needed for the chain of
    **      astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.
    **
    **      The various functions support different classes of observer and
    **      portions of the transformation chain:
    **
    **          functions         observer        transformation
    **
    **       iauApcg iauApcg13    geocentric      ICRS <-> GCRS
    **       iauApci iauApci13    terrestrial     ICRS <-> CIRS
    **       iauApco iauApco13    terrestrial     ICRS <-> observed
    **       iauApcs iauApcs13    space           ICRS <-> GCRS
    **       iauAper iauAper13    terrestrial     update Earth rotation
    **       iauApio iauApio13    terrestrial     CIRS <-> observed
    **
    **      Those with names ending in "13" use contemporary SOFA models to
    **      compute the various ephemerides.  The others accept ephemerides
    **      supplied by the caller.
    **
    **      The transformation from ICRS to GCRS covers space motion,
    **      parallax, light deflection, and aberration.  From GCRS to CIRS
    **      comprises frame bias and precession-nutation.  From CIRS to
    **      observed takes account of Earth rotation, polar motion, diurnal
    **      aberration and parallax (unless subsumed into the ICRS <-> GCRS
    **      transformation), and atmospheric refraction.
    **
    **  11) The context structure astrom produced by this function is used
    **      by iauAtioq and iauAtoiq.
    **
    **  Called:
    **     iauUtctai    UTC to TAI
    **     iauTaitt     TAI to TT
    **     iauUtcut1    UTC to UT1
    **     iauSp00      the TIO locator s', IERS 2000
    **     iauEra00     Earth rotation angle, IAU 2000
    **     iauRefco     refraction constants for given ambient conditions
    **     iauApio      astrometry parameters, CIRS-observed
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    astrom = ASTROM()
    s = _sofa.iauApio13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tc,rh,wl,astrom)
    if s != 0:
        warnings.warn(apio13_msg[s], UserWarning, 2)
    return astrom


# iauAtci13
try:
    _sofa.iauAtci13.argtypes = [
                            c_double, # rc /ICRS right ascension at J2000.0 (radians, Note 1)
                            c_double, # dc /ICRS declination at J2000.0 (radians, Note 1) 
                            c_double, # pr / RA proper motion (radians/year; Note 2)
                            c_double, # pd /Dec proper motion (radians/year)
                            c_double, # px /parallax (arcsec)
                            c_double, # rv /radial velocity (km/s, +ve if receding)
                            c_double, # date1 /TDB as a 2-part...
                            c_double, # date2 / ...Julian Date (Note 3)
                            POINTER(c_double), # ri /CIRS geocentric RA(radians)
                            POINTER(c_double), # di /CIRS geocentric Dec (radians)
                            POINTER(c_double)] # eo /equation of the origins (ERA-GST, Note 5)
except AttributeError:
    pass  
def atci13(rc,dc,pr,pd,px,rv,date1,date2):
    """
    /*
    **  - - - - - - - - - -
    **   i a u A t c i 1 3
    **  - - - - - - - - - -
    **
    **  Transform ICRS star data, epoch J2000.0, to CIRS.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     rc     double   ICRS right ascension at J2000.0 (radians, Note 1)
    **     dc     double   ICRS declination at J2000.0 (radians, Note 1)
    **     pr     double   RA proper motion (radians/year; Note 2)
    **     pd     double   Dec proper motion (radians/year)
    **     px     double   parallax (arcsec)
    **     rv     double   radial velocity (km/s, +ve if receding)
    **     date1  double   TDB as a 2-part...
    **     date2  double   ...Julian Date (Note 3)
    **
    **  Returned:
    **     ri,di  double*  CIRS geocentric RA,Dec (radians)
    **     eo     double*  equation of the origins (ERA-GST, Note 5)
    **
    **  Notes:
    **
    **  1) Star data for an epoch other than J2000.0 (for example from the
    **     Hipparcos catalog, which has an epoch of J1991.25) will require a
    **     preliminary call to iauPmsafe before use.
    **
    **  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
    **
    **  3) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  4) The available accuracy is better than 1 milliarcsecond, limited
    **     mainly by the precession-nutation model that is used, namely
    **     IAU 2000A/2006.  Very close to solar system bodies, additional
    **     errors of up to several milliarcseconds can occur because of
    **     unmodeled light deflection;  however, the Sun's contribution is
    **     taken into account, to first order.  The accuracy limitations of
    **     the SOFA function iauEpv00 (used to compute Earth position and
    **     velocity) can contribute aberration errors of up to
    **     5 microarcseconds.  Light deflection at the Sun's limb is
    **     uncertain at the 0.4 mas level.
    **
    **  5) Should the transformation to (equinox based) apparent place be
    **     required rather than (CIO based) intermediate place, subtract the
    **     equation of the origins from the returned right ascension:
    **     RA = RI - EO. (The iauAnp function can then be applied, as
    **     required, to keep the result in the conventional 0-2pi range.)
    **
    **  Called:
    **     iauApci13    astrometry parameters, ICRS-CIRS, 2013
    **     iauAtciq     quick ICRS to CIRS
    **
    **  This revision:   2017 March 12
    **
    **  SOFA release 2019-07-22
    **
    """
    ri = c_double()
    di = c_double()
    eo = c_double()
    _sofa.iauAtci13(rc,dc,pr,pd,px,rv,date1,date2,byref(ri),byref(di),byref(eo))
    return ri.value, di.value, eo.value


# iauAtciq
try:
    _sofa.iauAtciq.argtypes = [
                            c_double, # rc /ICRS right ascension at J2000.0 (radians)
                            c_double, # dc /ICRS declination at J2000.0 (radians) 
                            c_double, # pr / RA proper motion (radians/year; Note 3)
                            c_double, # pd /Dec proper motion (radians/year)
                            c_double, # px /parallax (arcsec)
                            c_double, # rv /radial velocity (km/s, +ve if receding)
                            POINTER(ASTROM), # astrom /star-independent astrometry parameters:
                            POINTER(c_double), # ri /CIRS geocentric RA(radians)
                            POINTER(c_double)] # di /CIRS geocentric Dec (radians))
except AttributeError:
    pass  
def atciq(rc,dc,pr,pd,px,rv,astrom):
    """
    /*
    **  - - - - - - - - -
    **   i a u A t c i q
    **  - - - - - - - - -
    **
    **  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
    **  star-independent astrometry parameters.
    **
    **  Use of this function is appropriate when efficiency is important and
    **  where many star positions are to be transformed for one date.  The
    **  star-independent parameters can be obtained by calling one of the
    **  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    **
    **  If the parallax and proper motions are zero the iauAtciqz function
    **  can be used instead.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     rc,dc  double     ICRS RA,Dec at J2000.0 (radians)
    **     pr     double     RA proper motion (radians/year; Note 3)
    **     pd     double     Dec proper motion (radians/year)
    **     px     double     parallax (arcsec)
    **     rv     double     radial velocity (km/s, +ve if receding)
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned:
    **     ri,di   double    CIRS RA,Dec (radians)
    **
    **  Notes:
    **
    **  1) All the vectors are with respect to BCRS axes.
    **
    **  2) Star data for an epoch other than J2000.0 (for example from the
    **     Hipparcos catalog, which has an epoch of J1991.25) will require a
    **     preliminary call to iauPmsafe before use.
    **
    **  3) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
    **
    **  Called:
    **     iauPmpx      proper motion and parallax
    **     iauLdsun     light deflection by the Sun
    **     iauAb        stellar aberration
    **     iauRxp       product of r-matrix and pv-vector
    **     iauC2s       p-vector to spherical
    **     iauAnp       normalize angle into range 0 to 2pi
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ri = c_double()
    di = c_double()
    _sofa.iauAtciq(rc,dc,pr,pd,px,rv,astrom,byref(ri),byref(di))
    return ri.value, di.value


#iauAtciqn
try:
    _sofa.iauAtciqn.argtypes = [
                            c_double, # rc /ICRS right ascension at J2000.0 (radians)
                            c_double, # dc /ICRS declination at J2000.0 (radians) 
                            c_double, # pr / RA proper motion (radians/year; Note 3)
                            c_double, # pd /Dec proper motion (radians/year)
                            c_double, # px /parallax (arcsec)
                            c_double, # rv /radial velocity (km/s, +ve if receding)
                            POINTER(ASTROM), # astrom /star-independent astrometry parameters:
                            c_int, # n /number of bodies (Note 3)
                            POINTER(LDBODY), # iauLDBODY[n] / data for each of the n bodies (Notes 3,4):
                            POINTER(c_double), # ri /CIRS geocentric RA(radians)
                            POINTER(c_double)] # di /CIRS geocentric Dec (radians))
except AttributeError:
    pass  
def atciqn(rc,dc,pr,pd,px,rv,astrom,n,ldbody):
    """
     Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed  star-independent astrometry parameters plus a list of light- deflecting bodies.

    Parameters
    ----------
    rc : double
        ICRS RA,Dec at J2000.0 (radians).
    dc : double
        ICRS RA,Dec at J2000.0 (radians).
    pr : double
        RA proper motion (radians/year; Note 3).
    pd : double
        Dec proper motion (radians/year).
    px : double
        parallax (arcsec).
    rv : double
        radial velocity (km/s, +ve if receding).
    astrom : TYPE
        star-independent astrometry parameters:.
    n : iauASTROM*
        number of bodies (Note 3).
    ldbody : iauLDBODY[n]
        data for each of the n bodies (Notes 3,4):.

    Returns
    -------
    ri : double
        CIRS RA,Dec (radians).
    di : double
        CIRS RA,Dec (radians).
/*
**  - - - - - - - - - -
**   i a u A t c i q n
**  - - - - - - - - - -
**
**  Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
**  star-independent astrometry parameters plus a list of light-
**  deflecting bodies.
**
**  Use of this function is appropriate when efficiency is important and
**  where many star positions are to be transformed for one date.  The
**  star-independent parameters can be obtained by calling one of the
**  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
**
**
**  If the only light-deflecting body to be taken into account is the
**  Sun, the iauAtciq function can be used instead.  If in addition the
**  parallax and proper motions are zero, the iauAtciqz function can be
**  used.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rc,dc  double       ICRS RA,Dec at J2000.0 (radians)
**     pr     double       RA proper motion (radians/year; Note 3)
**     pd     double       Dec proper motion (radians/year)
**     px     double       parallax (arcsec)
**     rv     double       radial velocity (km/s, +ve if receding)
**     astrom iauASTROM*   star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       longitude + s' (radians)
**      xpl    double       polar motion xp wrt local meridian (radians)
**      ypl    double       polar motion yp wrt local meridian (radians)
**      sphi   double       sine of geodetic latitude
**      cphi   double       cosine of geodetic latitude
**      diurab double       magnitude of diurnal aberration vector
**      eral   double       "local" Earth rotation angle (radians)
**      refa   double       refraction constant A (radians)
**      refb   double       refraction constant B (radians)
**      n     int           number of bodies (Note 3)
**      b     iauLDBODY[n] data for each of the n bodies (Notes 3,4):
**       bm    double        mass of the body (solar masses, Note 5)
**       dl    double        deflection limiter (Note 6)
**       pv    [2][3]        barycentric PV of the body (au, au/day)
**
**  Returned:
**     ri,di   double    CIRS RA,Dec (radians)
**
**  Notes:
**
**  1) Star data for an epoch other than J2000.0 (for example from the
**     Hipparcos catalog, which has an epoch of J1991.25) will require a
**     preliminary call to iauPmsafe before use.
**
**  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
**
**  3) The struct b contains n entries, one for each body to be
**     considered.  If n = 0, no gravitational light deflection will be
**     applied, not even for the Sun.
**
**  4) The struct b should include an entry for the Sun as well as for
**     any planet or other body to be taken into account.  The entries
**     should be in the order in which the light passes the body.
**
**  5) In the entry in the b struct for body i, the mass parameter
**     b[i].bm can, as required, be adjusted in order to allow for such
**     effects as quadrupole field.
**
**  6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
**     the angular separation (in radians) between star and body at
**     which limiting is applied.  As phi shrinks below the chosen
**     threshold, the deflection is artificially reduced, reaching zero
**     for phi = 0.   Example values suitable for a terrestrial
**     observer, together with masses, are as follows:
**
**        body i     b[i].bm        b[i].dl
**
**        Sun        1.0            6e-6
**        Jupiter    0.00095435     3e-9
**        Saturn     0.00028574     3e-10
**
**  7) For efficiency, validation of the contents of the b array is
**     omitted.  The supplied masses must be greater than zero, the
**     position and velocity vectors must be right, and the deflection
**     limiter greater than zero.
**
**  Called:
**     iauPmpx      proper motion and parallax
**     iauLdn       light deflection by n bodies
**     iauAb        stellar aberration
**     iauRxp       product of r-matrix and pv-vector
**     iauC2s       p-vector to spherical
**     iauAnp       normalize angle into range 0 to 2pi
**
**  This revision:   2013 October 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ri = c_double()
    di = c_double()
    _sofa.iauAtciqn(rc,dc,pr,pd,px,rv,astrom,n,ldbody,byref(ri),byref(di))
    return ri.value, di.value



# iauAtciqz
try:
    _sofa.iauAtciqz.argtypes = [
                            c_double, # rc /
                            c_double, # dc /ICRS declination at J2000.0 (radians) 
                            POINTER(ASTROM), # astrom /star-independent astrometry parameters:
                            POINTER(c_double), # ri /CIRS geocentric RA(radians)
                            POINTER(c_double) # di /CIRS geocentric Dec (radians))
                            ]
except AttributeError: 
    pass  
def atciqz(rc,dc,astrom):
    """
    

    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians).
    dc : float
        ICRS declination at J2000.0 (radians).
    astrom : structure/class
        CIRS geocentric Dec (radians)).

    Returns
    -------
    ri:float
        CIRS geocentric RA(radians)
    di:float
        CIRS geocentric Dec (radians)).

                                                      /*
    **  - - - - - - - - - -
    **   i a u A t c i q z
    **  -  - - - - - - - - -
    **
    **  Quick ICRS to CIRS transformation, given precomputed star-
    **  independent astrometry parameters, and assuming zero parallax and
    **  proper motion.
    **
    **  Use of this function is appropriate when efficiency is important and
    **  where many star positions are to be transformed for one date.  The
    **  star-independent parameters can be obtained by calling one of the
    **  functions iauApci[13], iauApcg[13], iauApco[13] or iauApcs[13].
    **
    **  The corresponding function for the case of non-zero parallax and
    **  proper motion is iauAtciq.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given: 
    **     rc,dc  double     ICRS astrometric RA,Dec (radians)
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned:
    **     ri,di  double     CIRS RA,Dec (radians)
    **
    **  Note:
    **
    **     All the vectors are with respect to BCRS axes.
    **
    **  References:
    **
    **     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    **     the Astronomical Almanac, 3rd ed., University Science Books
    **     (2013).
    **
    **     Klioner, Sergei A., "A practical relativistic model for micro-
    **     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
    **
    **  Called:
    **     iauS2c       spherical coordinates to unit vector
    **     iauLdsun     light deflection due to Sun
    **     iauAb        stellar aberration
    **     iauRxp       product of r-matrix and p-vector
    **     iauC2s       p-vector to spherical
    **     iauAnp       normalize angle into range +/- pi
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ri = c_double()
    di = c_double()
    _sofa.iauAtciqz(rc,dc,astrom,byref(ri),byref(di))
    return ri.value, di.value


 # iauAtco13
_sofa.iauAtco13.argtypes = [
    c_double,  # rc \RA(radians)
    c_double,  # dc \DEC(radians)
    c_double,  # pr \RA proper motion(radians/year)
    c_double,  # pd \DEC proper motion(radians/year)
    c_double,  # px \parallax(arcsec)
    c_double,  # rv \radial velocity(km/s,positive if receding)
    c_double,  # utc1 \the first part of UTC
    c_double,  # utc2 \the second part of UTC
    c_double,  # dut1 \UT1-UTC(seconds,)
    c_double,  # elong \longitude (raidans,east-positive)
    c_double,  # phi \geodetic latitude(radians)
    c_double,  # hm \height above ellipsoid(m)
    c_double,  # xp \X polar motion coordinates(radians)
    c_double,  # yp \Y polar motion coordinates(radians)
    c_double,  # phpa \ pressure at the observer(hPa=mB)
    c_double,  # tc \ambient temperature at the observer(deg C)
    c_double,  # rh \relative humidity at the observer(range 0-1)
    c_double,  # wl \wavelength(micrometers)
    POINTER(c_double),  # aob \ observed azimuth (radians:N=0degree,E=90degree)
    POINTER(c_double),  # zob \ observed zenith distance(radians)
    POINTER(c_double),  # hob \ observed hour angle(radians)
    POINTER(c_double),  # dob \.observed declination(radians)
    POINTER(c_double),  # rob \ observed right ascension(CIO-based,radians)
    POINTER(c_double),  # eo \ equation of the origins(ERA-GST)
]


def atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp,
           phpa, tc, rh, wl):
    """
    Parameters
    ----------
    rc : float
        ICRS right ascension at J2000.0 (radians, Note 1).
    dc : float
        ICRS right ascension at J2000.0 (radians, Note 1).
    pr : float
        RA proper motion (radians/year; Note 2)
    pd : float
        Dec proper motion (radians/year)
    px : float
         parallax (arcsec)
    rv : float
        radial velocity (km/s, +ve if receding)
    utc1 : float
        UTC as a 2-part...
    utc2 : float
        ...quasi Julian Date (Notes 3-4)
    dut1 : float
        UT1-UTC (seconds, Note 5)
    elong : float
        longitude (radians, east +ve, Note 6)
    phi : float
        latitude (geodetic, radians, Note 6)
    hm : float
        height above ellipsoid (m, geodetic, Notes 6,8)
    xp : float
       polar motion coordinates (radians, Note 7)
    yp : float
        polar motion coordinates (radians, Note 7)
    phpa : float
        pressure at the observer (hPa = mB, Note 8)
    tc : float
        ambient temperature at the observer (deg C)
    rh : float
        relative humidity at the observer (range 0-1)
    wl : float
        wavelength (micrometers, Note 9)

    Returns
    -------
    aob float
        observed azimuth (radians: N=0,E=90)
    zob : float
        observed zenith distance (radians)
    hob : float
        observed hour angle (radians)
    dob : float
        observed declination (radians)
    rob : float
        observed right ascension (CIO-based, radians)
    eo : float
        equation of the origins (ERA-GST)
        
        /*
    **  - - - - - - - - - -
    **   i a u A t c o 1 3
    **  - - - - - - - - - -
    **
    **  ICRS RA,Dec to observed place.  The caller supplies UTC, site
    **  coordinates, ambient air conditions and observing wavelength.
    **
    **  SOFA models are used for the Earth ephemeris, bias-precession-
    **  nutation, Earth orientation and refraction.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     rc,dc  double   ICRS right ascension at J2000.0 (radians, Note 1)
    **     pr     double   RA proper motion (radians/year; Note 2)
    **     pd     double   Dec proper motion (radians/year)
    **     px     double   parallax (arcsec)
    **     rv     double   radial velocity (km/s, +ve if receding)
    **     utc1   double   UTC as a 2-part...
    **     utc2   double   ...quasi Julian Date (Notes 3-4)
    **     dut1   double   UT1-UTC (seconds, Note 5)
    **     elong  double   longitude (radians, east +ve, Note 6)
    **     phi    double   latitude (geodetic, radians, Note 6)
    **     hm     double   height above ellipsoid (m, geodetic, Notes 6,8)
    **     xp,yp  double   polar motion coordinates (radians, Note 7)
    **     phpa   double   pressure at the observer (hPa = mB, Note 8)
    **     tc     double   ambient temperature at the observer (deg C)
    **     rh     double   relative humidity at the observer (range 0-1)
    **     wl     double   wavelength (micrometers, Note 9)
    **
    **  Returned:
    **     aob    double*  observed azimuth (radians: N=0,E=90)
    **     zob    double*  observed zenith distance (radians)
    **     hob    double*  observed hour angle (radians)
    **     dob    double*  observed declination (radians)
    **     rob    double*  observed right ascension (CIO-based, radians)
    **     eo     double*  equation of the origins (ERA-GST)
    **
    **  Returned (function value):
    **            int      status: +1 = dubious year (Note 4)
    **                              0 = OK
    **                             -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  Star data for an epoch other than J2000.0 (for example from the
    **      Hipparcos catalog, which has an epoch of J1991.25) will require
    **      a preliminary call to iauPmsafe before use.
    **
    **  2)  The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
    **
    **  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  4)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the
    **      future to be trusted.  See iauDat for further details.
    **
    **  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  6)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many
    **      applications, xp and yp can be set to zero.
    **
    **  8)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB),
    **      is available, an adequate estimate of hm can be obtained from
    **      the expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to
    **      the pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  9)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  10) The accuracy of the result is limited by the corrections for
    **      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **      Providing the meteorological parameters are known accurately and
    **      there are no gross local effects, the predicted observed
    **      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **      (radio) for a zenith distance of less than 70 degrees, better
    **      than 30 arcsec (optical or radio) at 85 degrees and better
    **      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **      Without refraction, the complementary functions iauAtco13 and
    **      iauAtoc13 are self-consistent to better than 1 microarcsecond
    **      all over the celestial sphere.  With refraction included,
    **      consistency falls off at high zenith distances, but is still
    **      better than 0.05 arcsec at 85 degrees.
    **
    **  11) "Observed" Az,ZD means the position that would be seen by a
    **      perfect geodetically aligned theodolite.  (Zenith distance is
    **      used rather than altitude in order to reflect the fact that no
    **      allowance is made for depression of the horizon.)  This is
    **      related to the observed HA,Dec via the standard rotation, using
    **      the geodetic latitude (corrected for polar motion), while the
    **      observed HA and RA are related simply through the Earth rotation
    **      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    **      means the position that would be seen by a perfect equatorial
    **      with its polar axis aligned to the Earth's axis of rotation.
    **
    **  12) It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  Called:
    **     iauApco13    astrometry parameters, ICRS-observed, 2013
    **     iauAtciq     quick ICRS to CIRS
    **     iauAtioq     quick CIRS to observed
    **
    **  This revision:   2016 February 2
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    aob = c_double()
    zob = c_double()
    hob = c_double()
    dob = c_double()
    rob = c_double()
    eo = c_double()
    _sofa.iauAtco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm,
                    xp, yp, phpa, tc, rh, wl, byref(aob), byref(zob), byref(hob), byref(dob), byref(rob), byref(eo))
    return aob.value, zob.value, hob.value, dob.value, rob.value, eo.value


# iauAtic13
try:
    _sofa.iauAtic13.argtypes = [
                            c_double, # ri
                            c_double, # di  
                            c_double, # date1
                            c_double, # date2
                            POINTER(c_double), # rc
                            POINTER(c_double), # dc
                            POINTER(c_double) # eo
                            ] 
except AttributeError:
    pass  
def atic13(ri,di,date1,date2):
    """
    Parameters
    ----------
    ri : float
        CIRS geocentric RA,Dec (radians).
    di : float
        CIRS geocentric RA,Dec (radians).
    date1 : TYPE
        TDB as a 2-part....
    date2 : TYPE
        ...Julian Date (Note 1).

    Returns
    -------
    rc : float
        ICRS astrometric RA,Dec (radians).
    dc : float
        ICRS astrometric RA,Dec (radians).
    eo : float
        equation of the origins (ERA-GST, Note 4).
    /*
    **  - - - - - - - - - -
    **   i a u A t i c 1 3
    **  - - - - - - - - - -
    **
    **  Transform star RA,Dec from geocentric CIRS to ICRS astrometric.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ri,di  double  CIRS geocentric RA,Dec (radians)
    **     date1  double  TDB as a 2-part...
    **     date2  double  ...Julian Date (Note 1)
    **
    **  Returned:
    **     rc,dc  double  ICRS astrometric RA,Dec (radians)
    **     eo     double  equation of the origins (ERA-GST, Note 4)
    **
    **  Notes:
    **
    **  1) The TDB date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TDB)=2450123.7 could be expressed in any of these ways, among
    **     others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.  For most
    **     applications of this function the choice will not be at all
    **     critical.
    **
    **     TT can be used instead of TDB without any significant impact on
    **     accuracy.
    **
    **  2) Iterative techniques are used for the aberration and light
    **     deflection corrections so that the functions iauAtic13 (or
    **     iauAticq) and iauAtci13 (or iauAtciq) are accurate inverses;
    **     even at the edge of the Sun's disk the discrepancy is only about
    **     1 nanoarcsecond.
    **
    **  3) The available accuracy is better than 1 milliarcsecond, limited
    **     mainly by the precession-nutation model that is used, namely
    **     IAU 2000A/2006.  Very close to solar system bodies, additional
    **     errors of up to several milliarcseconds can occur because of
    **     unmodeled light deflection;  however, the Sun's contribution is
    **     taken into account, to first order.  The accuracy limitations of
    **     the SOFA function iauEpv00 (used to compute Earth position and
    **     velocity) can contribute aberration errors of up to
    **     5 microarcseconds.  Light deflection at the Sun's limb is
    **     uncertain at the 0.4 mas level.
    **
    **  4) Should the transformation to (equinox based) J2000.0 mean place
    **     be required rather than (CIO based) ICRS coordinates, subtract the
    **     equation of the origins from the returned right ascension:
    **     RA = RI - EO.  (The iauAnp function can then be applied, as
    **     required, to keep the result in the conventional 0-2pi range.)
    **
    **  Called:
    **     iauApci13    astrometry parameters, ICRS-CIRS, 2013
    **     iauAticq     quick CIRS to ICRS astrometric
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    rc = c_double()
    dc = c_double()
    eo = c_double()
    _sofa.iauAtic13(ri,di,date1,date2,byref(rc),byref(dc),byref(eo))
    return rc.value, dc.value, eo.value


# iauAticq
try:
    _sofa.iauAticq.argtypes = [
                            c_double, # ri //CIRS RA at J2000.0 (radians)
                            c_double, # di /CIRS DEC at J2000.0 (radians) 
                            POINTER(ASTROM), # astrom /star-independent astrometry parameters:
                            POINTER(c_double), # rc /ICRS astrometric RA,Dec (radians)
                            POINTER(c_double) # dc /ICRS astrometric RA,Dec (radians)
                            ]
except AttributeError: 
    pass  
def aticq(ri,di,astrom):
    """
    Parameters
    ----------
    ri : float
        CIRS RA,Dec (radians).
    di : float
        CIRS RA,Dec (radians).
    astrom : structure
       star-independent astrometry parameters:

    Returns
    -------
    rc : float
        ICRS astrometric RA,Dec (radians).
    dc : float
        ICRS astrometric RA,Dec (radians).
    /*
    **  - - - - - - - - -
    **   i a u A t i c q
    **  - - - - - - - - -
    **
    **  Quick CIRS RA,Dec to ICRS astrometric place, given the star-
    **  independent astrometry parameters.
    **
    **  Use of this function is appropriate when efficiency is important and
    **  where many star positions are all to be transformed for one date.
    **  The star-independent astrometry parameters can be obtained by
    **  calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
    **  or iauApcs[13].
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ri,di  double     CIRS RA,Dec (radians)
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned:
    **     rc,dc  double     ICRS astrometric RA,Dec (radians)
    **
    **  Notes:
    **
    **  1) Only the Sun is taken into account in the light deflection
    **     correction.
    **
    **  2) Iterative techniques are used for the aberration and light
    **     deflection corrections so that the functions iauAtic13 (or
    **     iauAticq) and iauAtci13 (or iauAtciq) are accurate inverses;
    **     even at the edge of the Sun's disk the discrepancy is only about
    **     1 nanoarcsecond.
    **
    **  Called:
    **     iauS2c       spherical coordinates to unit vector
    **     iauTrxp      product of transpose of r-matrix and p-vector
    **     iauZp        zero p-vector
    **     iauAb        stellar aberration
    **     iauLdsun     light deflection by the Sun
    **     iauC2s       p-vector to spherical
    **     iauAnp       normalize angle into range +/- pi
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    rc = c_double()
    dc = c_double()
    _sofa.iauAticq(ri,di,astrom,byref(rc),byref(dc))
    return rc.value, dc.value

#iauAticqn
try:
    _sofa.iauAticqn.argtypes = [
                            c_double, # ri /ICRS right ascension at J2000.0 (radians)
                            c_double, # di /ICRS declination at J2000.0 (radians) 
                            POINTER(ASTROM), # astrom /star-independent astrometry parameters:
                            c_int, # n /number of bodies (Note 3)
                            POINTER(LDBODY), # iauLDBODY[n] / data for each of the n bodies (Notes 3,4):
                            POINTER(c_double), # rc /CIRS geocentric RA(radians)
                            POINTER(c_double)] # dc /CIRS geocentric Dec (radians))
except AttributeError:
    pass  
def aticqn(ri,di,astrom,n,ldbody):
    """
    Quick CIRS to ICRS astrometric place transformation, given the star-independent astrometry parameters plus a list of light-deflecting  bodies.

    Parameters
    ----------
    ri : double
        CIRS RA,Dec (radians).
    di : double
        CIRS RA,Dec (radians).
    astrom : iauASTROM
        star-independent astrometry parameters:.
    n : int
        number of bodies (Note 3).
    ldbody : iauLDBODY[n]
        data for each of the n bodies (Notes 3,4):.

    Returns
    -------
    rc : double
        ICRS astrometric RA,Dec (radians).
    dc : double
        ICRS astrometric RA,Dec (radians).
/*
**  - - - - - - - - - -
**   i a u A t i c q n
**  - - - - - - - - - -
**
**  Quick CIRS to ICRS astrometric place transformation, given the star-
**  independent astrometry parameters plus a list of light-deflecting
**  bodies.
**
**  Use of this function is appropriate when efficiency is important and
**  where many star positions are all to be transformed for one date.
**  The star-independent astrometry parameters can be obtained by
**  calling one of the functions iauApci[13], iauApcg[13], iauApco[13]
**  or iauApcs[13].
*
*  If the only light-deflecting body to be taken into account is the
*  Sun, the iauAticq function can be used instead.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     ri,di  double      CIRS RA,Dec (radians)
**     astrom iauASTROM*  star-independent astrometry parameters:
**      pmt    double       PM time interval (SSB, Julian years)
**      eb     double[3]    SSB to observer (vector, au)
**      eh     double[3]    Sun to observer (unit vector)
**      em     double       distance from Sun to observer (au)
**      v      double[3]    barycentric observer velocity (vector, c)
**      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
**      bpn    double[3][3] bias-precession-nutation matrix
**      along  double       longitude + s' (radians)
**      xpl    double       polar motion xp wrt local meridian (radians)
**      ypl    double       polar motion yp wrt local meridian (radians)
**      sphi   double       sine of geodetic latitude
**      cphi   double       cosine of geodetic latitude
**      diurab double       magnitude of diurnal aberration vector
**      eral   double       "local" Earth rotation angle (radians)
**      refa   double       refraction constant A (radians)
**      refb   double       refraction constant B (radians)
**      n     int           number of bodies (Note 3)
**      b     iauLDBODY[n] data for each of the n bodies (Notes 3,4):
**       bm    double       mass of the body (solar masses, Note 5)
**       dl    double       deflection limiter (Note 6)
**       pv    [2][3]       barycentric PV of the body (au, au/day)
**
**  Returned:
**     rc,dc  double     ICRS astrometric RA,Dec (radians)
**
**  Notes:
**
**  1) Iterative techniques are used for the aberration and light
**     deflection corrections so that the functions iauAticqn and
**     iauAtciqn are accurate inverses; even at the edge of the Sun's
**     disk the discrepancy is only about 1 nanoarcsecond.
**
**  2) If the only light-deflecting body to be taken into account is the
**     Sun, the iauAticq function can be used instead.
**
**  3) The struct b contains n entries, one for each body to be
**     considered.  If n = 0, no gravitational light deflection will be
**     applied, not even for the Sun.
**
**  4) The struct b should include an entry for the Sun as well as for
**     any planet or other body to be taken into account.  The entries
**     should be in the order in which the light passes the body.
**
**  5) In the entry in the b struct for body i, the mass parameter
**     b[i].bm can, as required, be adjusted in order to allow for such
**     effects as quadrupole field.
**
**  6) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
**     the angular separation (in radians) between star and body at
**     which limiting is applied.  As phi shrinks below the chosen
**     threshold, the deflection is artificially reduced, reaching zero
**     for phi = 0.   Example values suitable for a terrestrial
**     observer, together with masses, are as follows:
**
**        body i     b[i].bm        b[i].dl
**
**        Sun        1.0            6e-6
**        Jupiter    0.00095435     3e-9
**        Saturn     0.00028574     3e-10
**
**  7) For efficiency, validation of the contents of the b array is
**     omitted.  The supplied masses must be greater than zero, the
**     position and velocity vectors must be right, and the deflection
**     limiter greater than zero.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauZp        zero p-vector
**     iauAb        stellar aberration
**     iauLdn       light deflection by n bodies
**     iauC2s       p-vector to spherical
**     iauAnp       normalize angle into range +/- pi
**
**  This revision:   2018 December 13
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc = c_double()
    dc = c_double()
    _sofa.iauAticqn(ri,di,astrom,n,ldbody,byref(rc),byref(dc))
    return rc.value, dc.value

 # iauAtio13
atio13_msg = {0: 'OK', # Unused
                +1:'dubious year (Note 2)',
                -1:'unacceptable date',
                }
_sofa.iauAtio13.argtypes = [
    c_double,  # ri \RA(radians)
    c_double,  # di \DEC(radians)
    c_double,  # utc1 \the first part of UTC
    c_double,  # utc2 \the second part of UTC
    c_double,  # dut1 \UT1-UTC(seconds,)
    c_double,  # elong \longitude (raidans,east-positive)
    c_double,  # phi \geodetic latitude(radians)
    c_double,  # hm \height above ellipsoid(m)
    c_double,  # xp \X polar motion coordinates(radians)
    c_double,  # yp \Y polar motion coordinates(radians)
    c_double,  # phpa \ pressure at the observer(hPa=mB)
    c_double,  # tc \ambient temperature at the observer(deg C)
    c_double,  # rh \relative humidity at the observer(range 0-1)
    c_double,  # wl \wavelength(micrometers)
    POINTER(c_double),  # aob \ observed azimuth (radians:N=0degree,E=90degree)
    POINTER(c_double),  # zob \ observed zenith distance(radians)
    POINTER(c_double),  # hob \ observed hour angle(radians)
    POINTER(c_double),  # dob \.observed declination(radians)
    POINTER(c_double),  # rob \ observed right ascension(CIO-based,radians)
]
def atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl):
    """
    Parameters
    ----------
    ri : float
        CIRS right ascension (CIO-based, radians).
    di : float
        CIRS declination (radians).
    utc1 : float
        UTC as a 2-part...
    utc2 : float
        ...quasi Julian Date (Notes 1,2)
    dut1 : float
       UT1-UTC (seconds, Note 3)
    elong : float
        longitude (radians, east +ve, Note 4)
    phi : float
        geodetic latitude (radians, Note 4)
    hm : float
        height above ellipsoid (m, geodetic Notes 4,6)
    xp : float
        polar motion coordinates (radians, Note 5)
    yp : float
        polar motion coordinates (radians, Note 5).
    phpa : float
        pressure at the observer (hPa = mB, Note 6).
    tc : float
        ambient temperature at the observer (deg C).
    rh : float
        relative humidity at the observer (range 0-1).
    wl : float
        wavelength (micrometers, Note 7).

    Returns
    -------
    aob float
        observed azimuth (radians: N=0,E=90)
    zob : float
        observed zenith distance (radians)
    hob : float
        observed hour angle (radians)
    dob : float
        observed declination (radians)
    rob : float
        observed right ascension (CIO-based, radians)
    eo : float
        equation of the origins (ERA-GST)
    /*
    **  - - - - - - - - - -
    **   i a u A t i o 1 3
    **  - - - - - - - - - -
    **
    **  CIRS RA,Dec to observed place.  The caller supplies UTC, site
    **  coordinates, ambient air conditions and observing wavelength.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ri     double   CIRS right ascension (CIO-based, radians)
    **     di     double   CIRS declination (radians)
    **     utc1   double   UTC as a 2-part...
    **     utc2   double   ...quasi Julian Date (Notes 1,2)
    **     dut1   double   UT1-UTC (seconds, Note 3)
    **     elong  double   longitude (radians, east +ve, Note 4)
    **     phi    double   geodetic latitude (radians, Note 4)
    **     hm     double   height above ellipsoid (m, geodetic Notes 4,6)
    **     xp,yp  double   polar motion coordinates (radians, Note 5)
    **     phpa   double   pressure at the observer (hPa = mB, Note 6)
    **     tc     double   ambient temperature at the observer (deg C)
    **     rh     double   relative humidity at the observer (range 0-1)
    **     wl     double   wavelength (micrometers, Note 7)
    **
    **  Returned:
    **     aob    double*  observed azimuth (radians: N=0,E=90)
    **     zob    double*  observed zenith distance (radians)
    **     hob    double*  observed hour angle (radians)
    **     dob    double*  observed declination (radians)
    **     rob    double*  observed right ascension (CIO-based, radians)
    **
    **  Returned (function value):
    **            int      status: +1 = dubious year (Note 2)
    **                              0 = OK
    **                             -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  2)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the
    **      future to be trusted.  See iauDat for further details.
    **
    **  3)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  4)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  5)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many
    **      applications, xp and yp can be set to zero.
    **
    **  6)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB), is
    **      available, an adequate estimate of hm can be obtained from the
    **      expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to
    **      the pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  7)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  8)  "Observed" Az,ZD means the position that would be seen by a
    **      perfect geodetically aligned theodolite.  (Zenith distance is
    **      used rather than altitude in order to reflect the fact that no
    **      allowance is made for depression of the horizon.)  This is
    **      related to the observed HA,Dec via the standard rotation, using
    **      the geodetic latitude (corrected for polar motion), while the
    **      observed HA and RA are related simply through the Earth rotation
    **      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    **      means the position that would be seen by a perfect equatorial
    **      with its polar axis aligned to the Earth's axis of rotation.
    **
    **  9)  The accuracy of the result is limited by the corrections for
    **      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **      Providing the meteorological parameters are known accurately and
    **      there are no gross local effects, the predicted astrometric
    **      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **      (radio) for a zenith distance of less than 70 degrees, better
    **      than 30 arcsec (optical or radio) at 85 degrees and better
    **      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **  10) The complementary functions iauAtio13 and iauAtoi13 are self-
    **      consistent to better than 1 microarcsecond all over the
    **      celestial sphere.
    **
    **  11) It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  Called:
    **     iauApio13    astrometry parameters, CIRS-observed, 2013
    **     iauAtioq     quick CIRS to observed
    **
    **  This revision:   2016 February 2
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    aob = c_double()
    zob = c_double()
    hob = c_double()
    dob = c_double()
    rob = c_double()
    s=_sofa.iauAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, 
                      phpa, tc, rh, wl,byref(aob),byref(zob),byref(hob),byref(dob),byref(rob))
    if s != 0:
        warnings.warn(atio13_msg[s], UserWarning, 2)
    return aob.value, zob.value, hob.value, dob.value, rob.value

# iauAtioq
_sofa.iauAtioq.argtypes = [
    c_double,  # ri \RA(radians)
    c_double,  # di \DEC(radians)
    POINTER(ASTROM),  # astrom
    POINTER(c_double),  # aob \ observed azimuth (radians:N=0degree,E=90degree)
    POINTER(c_double),  # zob \ observed zenith distance(radians)
    POINTER(c_double),  # hob \ observed hour angle(radians)
    POINTER(c_double),  # dob \.observed declination(radians)
    POINTER(c_double),  # rob \ observed right ascension(CIO-based,radians)
]
def atioq(ri, di, astrom):
    """
    

    Parameters
    ----------
    ri : TYPE
        CIRS right ascension.
    di : TYPE
        CIRS declination.
    astrom : TYPE
        star-independent astrometry parameters:.

    Returns
    -------
    **     aob    double*    observed azimuth (radians: N=0,E=90)
    **     zob    double*    observed zenith distance (radians)
    **     hob    double*    observed hour angle (radians)
    **     dob    double*    observed declination (radians)
    **     rob    double*    observed right ascension (CIO-based, radians)
    /*
    **  - - - - - - - - -
    **   i a u A t i o q
    **  - - - - - - - - -
    **
    **  Quick CIRS to observed place transformation.
    **
    **  Use of this function is appropriate when efficiency is important and
    **  where many star positions are all to be transformed for one date.
    **  The star-independent astrometry parameters can be obtained by
    **  calling iauApio[13] or iauApco[13].
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ri     double     CIRS right ascension
    **     di     double     CIRS declination
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned:
    **     aob    double*    observed azimuth (radians: N=0,E=90)
    **     zob    double*    observed zenith distance (radians)
    **     hob    double*    observed hour angle (radians)
    **     dob    double*    observed declination (radians)
    **     rob    double*    observed right ascension (CIO-based, radians)
    **
    **  Notes:
    **
    **  1) This function returns zenith distance rather than altitude in
    **     order to reflect the fact that no allowance is made for
    **     depression of the horizon.
    **
    **  2) The accuracy of the result is limited by the corrections for
    **     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **     Providing the meteorological parameters are known accurately and
    **     there are no gross local effects, the predicted observed
    **     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **     (radio) for a zenith distance of less than 70 degrees, better
    **     than 30 arcsec (optical or radio) at 85 degrees and better
    **     than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **     Without refraction, the complementary functions iauAtioq and
    **     iauAtoiq are self-consistent to better than 1 microarcsecond all
    **     over the celestial sphere.  With refraction included, consistency
    **     falls off at high zenith distances, but is still better than
    **     0.05 arcsec at 85 degrees.
    **
    **  3) It is advisable to take great care with units, as even unlikely
    **     values of the input parameters are accepted and processed in
    **     accordance with the models used.
    **
    **  4) The CIRS RA,Dec is obtained from a star catalog mean place by
    **     allowing for space motion, parallax, the Sun's gravitational lens
    **     effect, annual aberration and precession-nutation.  For star
    **     positions in the ICRS, these effects can be applied by means of
    **     the iauAtci13 (etc.) functions.  Starting from classical "mean
    **     place" systems, additional transformations will be needed first.
    **
    **  5) "Observed" Az,El means the position that would be seen by a
    **     perfect geodetically aligned theodolite.  This is obtained from
    **     the CIRS RA,Dec by allowing for Earth orientation and diurnal
    **     aberration, rotating from equator to horizon coordinates, and
    **     then adjusting for refraction.  The HA,Dec is obtained by
    **     rotating back into equatorial coordinates, and is the position
    **     that would be seen by a perfect equatorial with its polar axis
    **     aligned to the Earth's axis of rotation.  Finally, the RA is
    **     obtained by subtracting the HA from the local ERA.
    **
    **  6) The star-independent CIRS-to-observed-place parameters in ASTROM
    **     may be computed with iauApio[13] or iauApco[13].  If nothing has
    **     changed significantly except the time, iauAper[13] may be used to
    **     perform the requisite adjustment to the astrom structure.
    **
    **  Called:
    **     iauS2c       spherical coordinates to unit vector
    **     iauC2s       p-vector to spherical
    **     iauAnp       normalize angle into range 0 to 2pi
    **
    **  This revision:   2016 March 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    aob = c_double()
    zob = c_double()
    hob = c_double()
    dob = c_double()
    rob = c_double()
    _sofa.iauAtioq(ri, di,astrom,byref(aob),byref(zob),byref(hob),byref(dob),byref(rob))
    return aob.value, zob.value, hob.value, dob.value, rob.value

# iauAtoc13
atoc13_msg = {0: 'OK', # Unused
                +1:'dubious year (Note 2)',
                -1:'unacceptable date',
                }
_sofa.iauAtoc13.argtypes = [
    POINTER(c_char),   # type, char[] 
    c_longdouble, # ob1
    c_longdouble, # ob2
    c_longdouble,  # utc1 \the first part of UTC
    c_longdouble,  # utc2 \the second part of UTC
    c_longdouble,  # dut1 \UT1-UTC(seconds,)
    c_longdouble,  # elong \longitude (raidans,east-positive)
    c_longdouble,  # phi \geodetic latitude(radians)
    c_longdouble,  # hm \height above ellipsoid(m)
    c_longdouble,  # xp \X polar motion coordinates(radians)
    c_longdouble,  # yp \Y polar motion coordinates(radians)
    c_longdouble,  # phpa \ pressure at the observer(hPa=mB)
    c_longdouble,  # tc \ambient temperature at the observer(deg C)
    c_longdouble,  # rh \relative humidity at the observer(range 0-1)
    c_longdouble,  # wl \wavelength(micrometers)
    POINTER(c_longdouble),  # rc \ ICRS astrometric RA,Dec (radians)
    POINTER(c_longdouble)  # dc \ ICRS astrometric RA,Dec (radians)
]
def atoc13(strtype, ob1,ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, 
           phpa, tc, rh, wl):
    """
    

    Parameters
    ----------
    strtype : str
        type of coordinates - "R", "H" or "A" (Notes 1,2)
    ob1 : TYPE
        observed Az, HA or RA (radians; Az is N=0,E=90)
    ob2 : TYPE
        observed ZD or Dec (radians)
    utc1 : TYPE
        UTC as a 2-part...
    utc2 : TYPE
        ...quasi Julian Date (Notes 3,4)
    dut1 : TYPE
        UT1-UTC (seconds, Note 5)
    elong : TYPE
        longitude (radians, east +ve, Note 6)
    phi : TYPE
        geodetic latitude (radians, Note 6)
    hm : TYPE
        height above ellipsoid (m, geodetic Notes 6,8).
    xp : TYPE
        polar motion coordinates (radians, Note 7
    yp : TYPE
        polar motion coordinates (radians, Note 7
    phpa : TYPE
        pressure at the observer (hPa = mB, Note 8).
    tc : TYPE
        ambient temperature at the observer (deg C).
    rh : TYPE
       relative humidity at the observer (range 0-1)
    wl : TYPE
       wavelength (micrometers, Note 9).

    Returns
    -------
    rc :TYPE
        ICRS astrometric RA (radians)
    dc: TYPE
        ICRS astrometric Dec (radians)
    /*
    **  - - - - - - - - - -
    **   i a u A t o c 1 3
    **  - - - - - - - - - -
    **
    **  Observed place at a groundbased site to to ICRS astrometric RA,Dec.
    **  The caller supplies UTC, site coordinates, ambient air conditions
    **  and observing wavelength.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     type   char[]   type of coordinates - "R", "H" or "A" (Notes 1,2)
    **     ob1    double   observed Az, HA or RA (radians; Az is N=0,E=90)
    **     ob2    double   observed ZD or Dec (radians)
    **     utc1   double   UTC as a 2-part...
    **     utc2   double   ...quasi Julian Date (Notes 3,4)
    **     dut1   double   UT1-UTC (seconds, Note 5)
    **     elong  double   longitude (radians, east +ve, Note 6)
    **     phi    double   geodetic latitude (radians, Note 6)
    **     hm     double   height above ellipsoid (m, geodetic Notes 6,8)
    **     xp,yp  double   polar motion coordinates (radians, Note 7)
    **     phpa   double   pressure at the observer (hPa = mB, Note 8)
    **     tc     double   ambient temperature at the observer (deg C)
    **     rh     double   relative humidity at the observer (range 0-1)
    **     wl     double   wavelength (micrometers, Note 9)
    **
    **  Returned:
    **     rc,dc  double   ICRS astrometric RA,Dec (radians)
    **
    **  Returned (function value):
    **            int      status: +1 = dubious year (Note 4)
    **                              0 = OK
    **                             -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  "Observed" Az,ZD means the position that would be seen by a
    **      perfect geodetically aligned theodolite.  (Zenith distance is
    **      used rather than altitude in order to reflect the fact that no
    **      allowance is made for depression of the horizon.)  This is
    **      related to the observed HA,Dec via the standard rotation, using
    **      the geodetic latitude (corrected for polar motion), while the
    **      observed HA and RA are related simply through the Earth rotation
    **      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    **      means the position that would be seen by a perfect equatorial
    **      with its polar axis aligned to the Earth's axis of rotation.
    **
    **  2)  Only the first character of the type argument is significant.
    **      "R" or "r" indicates that ob1 and ob2 are the observed right
    **      ascension and declination;  "H" or "h" indicates that they are
    **      hour angle (west +ve) and declination;  anything else ("A" or
    **      "a" is recommended) indicates that ob1 and ob2 are azimuth
    **      (north zero, east 90 deg) and zenith distance.
    **
    **  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  4)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the
    **      future to be trusted.  See iauDat for further details.
    **
    **  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  6)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many
    **      applications, xp and yp can be set to zero.
    **
    **  8)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB), is
    **      available, an adequate estimate of hm can be obtained from the
    **      expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to
    **      the pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  9)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  10) The accuracy of the result is limited by the corrections for
    **      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **      Providing the meteorological parameters are known accurately and
    **      there are no gross local effects, the predicted astrometric
    **      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **      (radio) for a zenith distance of less than 70 degrees, better
    **      than 30 arcsec (optical or radio) at 85 degrees and better
    **      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **      Without refraction, the complementary functions iauAtco13 and
    **      iauAtoc13 are self-consistent to better than 1 microarcsecond
    **      all over the celestial sphere.  With refraction included,
    **      consistency falls off at high zenith distances, but is still
    **      better than 0.05 arcsec at 85 degrees.
    **
    **  11) It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  Called:
    **     iauApco13    astrometry parameters, ICRS-observed
    **     iauAtoiq     quick observed to CIRS
    **     iauAticq     quick CIRS to ICRS
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    rc = c_longdouble()
    dc = c_longdouble()
    strr = strtype.encode()
    s = _sofa.iauAtoc13(strr,ob1,ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,byref(rc),byref(dc))
    if s != 0:
        warnings.warn(atoc13_msg[s], UserWarning, 2)
    return rc.value, dc.value


# iauAtoi13
atoi13_msg = {0: 'OK', # Unused
                +1:'dubious year (Note 2)',
                -1:'unacceptable date',
                }
_sofa.iauAtoi13.argtypes = [
    POINTER(c_char),   # type, char[] 
    c_longdouble, # ob1
    c_longdouble, # ob2
    c_longdouble,  # utc1 \the first part of UTC
    c_longdouble,  # utc2 \the second part of UTC
    c_longdouble,  # dut1 \UT1-UTC(seconds,)
    c_longdouble,  # elong \longitude (raidans,east-positive)
    c_longdouble,  # phi \geodetic latitude(radians)
    c_longdouble,  # hm \height above ellipsoid(m)
    c_longdouble,  # xp \X polar motion coordinates(radians)
    c_longdouble,  # yp \Y polar motion coordinates(radians)
    c_longdouble,  # phpa \ pressure at the observer(hPa=mB)
    c_longdouble,  # tc \ambient temperature at the observer(deg C)
    c_longdouble,  # rh \relative humidity at the observer(range 0-1)
    c_longdouble,  # wl \wavelength(micrometers)
    POINTER(c_longdouble),  # ri \ CIRS right ascension (CIO-based, radians)
    POINTER(c_longdouble)  # di \ CIRS declination (radians)
]
def atoi13(strtype, ob1,ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, 
           phpa, tc, rh, wl):
    """
        Parameters
    ----------
    strtype : str
        type of coordinates - "R", "H" or "A" (Notes 1,2)
    ob1 : TYPE
        observed Az, HA or RA (radians; Az is N=0,E=90)
    ob2 : TYPE
        observed ZD or Dec (radians)
    utc1 : TYPE
        UTC as a 2-part...
    utc2 : TYPE
        ...quasi Julian Date (Notes 3,4)
    dut1 : TYPE
        UT1-UTC (seconds, Note 5)
    elong : TYPE
        longitude (radians, east +ve, Note 6)
    phi : TYPE
        geodetic latitude (radians, Note 6)
    hm : TYPE
        height above ellipsoid (m, geodetic Notes 6,8).
    xp : TYPE
        polar motion coordinates (radians, Note 7
    yp : TYPE
        polar motion coordinates (radians, Note 7
    phpa : TYPE
        pressure at the observer (hPa = mB, Note 8).
    tc : TYPE
        ambient temperature at the observer (deg C).
    rh : TYPE
       relative humidity at the observer (range 0-1)
    wl : TYPE
       wavelength (micrometers, Note 9).

    Returns
    -------
    ri :TYPE
        CIRS right ascension (CIO-based, radians)
    di: TYPE
        CIRS declination (radians)

    Parameters
    ----------

    /*
    **  - - - - - - - - - -
    **   i a u A t o i 1 3
    **  - - - - - - - - - -
    **
    **  Observed place to CIRS.  The caller supplies UTC, site coordinates,
    **  ambient air conditions and observing wavelength.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     type   char[]   type of coordinates - "R", "H" or "A" (Notes 1,2)
    **     ob1    double   observed Az, HA or RA (radians; Az is N=0,E=90)
    **     ob2    double   observed ZD or Dec (radians)
    **     utc1   double   UTC as a 2-part...
    **     utc2   double   ...quasi Julian Date (Notes 3,4)
    **     dut1   double   UT1-UTC (seconds, Note 5)
    **     elong  double   longitude (radians, east +ve, Note 6)
    **     phi    double   geodetic latitude (radians, Note 6)
    **     hm     double   height above the ellipsoid (meters, Notes 6,8)
    **     xp,yp  double   polar motion coordinates (radians, Note 7)
    **     phpa   double   pressure at the observer (hPa = mB, Note 8)
    **     tc     double   ambient temperature at the observer (deg C)
    **     rh     double   relative humidity at the observer (range 0-1)
    **     wl     double   wavelength (micrometers, Note 9)
    **
    **  Returned:
    **     ri     double*  CIRS right ascension (CIO-based, radians)
    **     di     double*  CIRS declination (radians)
    **
    **  Returned (function value):
    **            int      status: +1 = dubious year (Note 2)
    **                              0 = OK
    **                             -1 = unacceptable date
    **
    **  Notes:
    **
    **  1)  "Observed" Az,ZD means the position that would be seen by a
    **      perfect geodetically aligned theodolite.  (Zenith distance is
    **      used rather than altitude in order to reflect the fact that no
    **      allowance is made for depression of the horizon.)  This is
    **      related to the observed HA,Dec via the standard rotation, using
    **      the geodetic latitude (corrected for polar motion), while the
    **      observed HA and RA are related simply through the Earth rotation
    **      angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    **      means the position that would be seen by a perfect equatorial
    **      with its polar axis aligned to the Earth's axis of rotation.
    **
    **  2)  Only the first character of the type argument is significant.
    **      "R" or "r" indicates that ob1 and ob2 are the observed right
    **      ascension and declination;  "H" or "h" indicates that they are
    **      hour angle (west +ve) and declination;  anything else ("A" or
    **      "a" is recommended) indicates that ob1 and ob2 are azimuth
    **      (north zero, east 90 deg) and zenith distance.
    **
    **  3)  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    **      convenient way between the two arguments, for example where utc1
    **      is the Julian Day Number and utc2 is the fraction of a day.
    **
    **      However, JD cannot unambiguously represent UTC during a leap
    **      second unless special measures are taken.  The convention in the
    **      present function is that the JD day represents UTC days whether
    **      the length is 86399, 86400 or 86401 SI seconds.
    **
    **      Applications should use the function iauDtf2d to convert from
    **      calendar date and time of day into 2-part quasi Julian Date, as
    **      it implements the leap-second-ambiguity convention just
    **      described.
    **
    **  4)  The warning status "dubious year" flags UTCs that predate the
    **      introduction of the time scale or that are too far in the
    **      future to be trusted.  See iauDat for further details.
    **
    **  5)  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    **      one second at the end of each positive UTC leap second,
    **      introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    **      practice is under review, and in the future UT1-UTC may grow
    **      essentially without limit.
    **
    **  6)  The geographical coordinates are with respect to the WGS84
    **      reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    **      longitude required by the present function is east-positive
    **      (i.e. right-handed), in accordance with geographical convention.
    **
    **  7)  The polar motion xp,yp can be obtained from IERS bulletins.  The
    **      values are the coordinates (in radians) of the Celestial
    **      Intermediate Pole with respect to the International Terrestrial
    **      Reference System (see IERS Conventions 2003), measured along the
    **      meridians 0 and 90 deg west respectively.  For many
    **      applications, xp and yp can be set to zero.
    **
    **  8)  If hm, the height above the ellipsoid of the observing station
    **      in meters, is not known but phpa, the pressure in hPa (=mB), is
    **      available, an adequate estimate of hm can be obtained from the
    **      expression
    **
    **            hm = -29.3 * tsl * log ( phpa / 1013.25 );
    **
    **      where tsl is the approximate sea-level air temperature in K
    **      (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    **      52).  Similarly, if the pressure phpa is not known, it can be
    **      estimated from the height of the observing station, hm, as
    **      follows:
    **
    **            phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    **
    **      Note, however, that the refraction is nearly proportional to
    **      the pressure and that an accurate phpa value is important for
    **      precise work.
    **
    **  9)  The argument wl specifies the observing wavelength in
    **      micrometers.  The transition from optical to radio is assumed to
    **      occur at 100 micrometers (about 3000 GHz).
    **
    **  10) The accuracy of the result is limited by the corrections for
    **      refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **      Providing the meteorological parameters are known accurately and
    **      there are no gross local effects, the predicted astrometric
    **      coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **      (radio) for a zenith distance of less than 70 degrees, better
    **      than 30 arcsec (optical or radio) at 85 degrees and better
    **      than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **      Without refraction, the complementary functions iauAtio13 and
    **      iauAtoi13 are self-consistent to better than 1 microarcsecond
    **      all over the celestial sphere.  With refraction included,
    **      consistency falls off at high zenith distances, but is still
    **      better than 0.05 arcsec at 85 degrees.
    **
    **  12) It is advisable to take great care with units, as even unlikely
    **      values of the input parameters are accepted and processed in
    **      accordance with the models used.
    **
    **  Called:
    **     iauApio13    astrometry parameters, CIRS-observed, 2013
    **     iauAtoiq     quick observed to CIRS
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ri = c_longdouble()
    di = c_longdouble()
    strr = strtype.encode()
    s = _sofa.iauAtoi13(strr,ob1,ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl,byref(ri),byref(di))
    if s != 0:
        warnings.warn(atoi13_msg[s], UserWarning, 2)
    return ri.value, di.value

#iauLdn
_sofa.iauLdn.argtypes = [
    c_int,   # n,
    POINTER(LDBODY),
    ndpointer(shape=3,dtype=float), # ob[3]
    ndpointer(shape=3,dtype=float), # sc[3]
    ndpointer(shape=3,dtype=float), # sn[3]
]
def ldn(n,ldbody,ob,sc):
    """
    For a star, apply light deflection by multiple solar-system bodies,  as part of transforming coordinate direction into natural direction.

    Parameters
    ----------
    n : int
        number of bodies (note 1).
    ldbody : iauLDBODY[n]
        data for each of the n bodies (Notes 1,2):.
    ob : double[3]
        barycentric position of the observer (au).
    sc : double[3]
        observer to star coord direction (unit vector).

    Returns
    -------
    sn : double[3]
        observer to deflected star (unit vector).
/*+
**  - - - - - - -
**   i a u L d n
**  - - - - - - -
**
**  For a star, apply light deflection by multiple solar-system bodies,
**  as part of transforming coordinate direction into natural direction.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     n    int           number of bodies (note 1)
**     b    iauLDBODY[n]  data for each of the n bodies (Notes 1,2):
**      bm   double         mass of the body (solar masses, Note 3)
**      dl   double         deflection limiter (Note 4)
**      pv   [2][3]         barycentric PV of the body (au, au/day)
**     ob   double[3]     barycentric position of the observer (au)
**     sc   double[3]     observer to star coord direction (unit vector)
**
**  Returned:
**     sn    double[3]      observer to deflected star (unit vector)
**
**  1) The array b contains n entries, one for each body to be
**     considered.  If n = 0, no gravitational light deflection will be
**     applied, not even for the Sun.
**
**  2) The array b should include an entry for the Sun as well as for
**     any planet or other body to be taken into account.  The entries
**     should be in the order in which the light passes the body.
**
**  3) In the entry in the b array for body i, the mass parameter
**     b[i].bm can, as required, be adjusted in order to allow for such
**     effects as quadrupole field.
**
**  4) The deflection limiter parameter b[i].dl is phi^2/2, where phi is
**     the angular separation (in radians) between star and body at
**     which limiting is applied.  As phi shrinks below the chosen
**     threshold, the deflection is artificially reduced, reaching zero
**     for phi = 0.   Example values suitable for a terrestrial
**     observer, together with masses, are as follows:
**
**        body i     b[i].bm        b[i].dl
**
**        Sun        1.0            6e-6
**        Jupiter    0.00095435     3e-9
**        Saturn     0.00028574     3e-10
**
**  5) For cases where the starlight passes the body before reaching the
**     observer, the body is placed back along its barycentric track by
**     the light time from that point to the observer.  For cases where
**     the body is "behind" the observer no such shift is applied.  If
**     a different treatment is preferred, the user has the option of
**     instead using the iauLd function.  Similarly, iauLd can be used
**     for cases where the source is nearby, not a star.
**
**  6) The returned vector sn is not normalized, but the consequential
**     departure from unit magnitude is always negligible.
**
**  7) The arguments sc and sn can be the same array.
**
**  8) For efficiency, validation is omitted.  The supplied masses must
**     be greater than zero, the position and velocity vectors must be
**     right, and the deflection limiter greater than zero.
**
**  Reference:
**
**     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
**     the Astronomical Almanac, 3rd ed., University Science Books
**     (2013), Section 7.2.4.
**
**  Called:
**     iauCp        copy p-vector
**     iauPdp       scalar product of two p-vectors
**     iauPmp       p-vector minus p-vector
**     iauPpsp      p-vector plus scaled p-vector
**     iauPn        decompose p-vector into modulus and direction
**     iauLd        light deflection by a solar-system body
**
**  This revision:   2017 March 16
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    sn = zeros(shape=3,dtype=float)
    _sofa.iauLdn(n,ldbody,ob,sc,sn)
    return sn

#iauAtoiq
_sofa.iauAtoiq.argtypes = [
    POINTER(c_char),   # strtype, char[] 
    c_double, # ob1
    c_double, # ob2
    POINTER(ASTROM),
    POINTER(c_longdouble),  # ri \ CIRS right ascension (CIO-based, radians)
    POINTER(c_longdouble)  # di \ CIRS declination (radians)
]
def atoiq(strtype, ob1,ob2, astrom):
    """
    

    Parameters
    ----------
    strtype : char
        type of coordinates: "R", "H" or "A" (Note 1)
    ob1 : TYPE
       observed Az, HA or RA (radians; Az is N=0,E=90)
    ob2 : TYPE
        observed ZD or Dec (radians)
    astrom : TYPE
        star-independent astrometry parameters:

    Returns
    -------
    ri :TYPE
        CIRS right ascension (CIO-based, radians)
    di: TYPE
        CIRS declination (radians)

    Parameters
    ----------
    /*
    **  - - - - - - - - -
    **   i a u A t o i q
    **  - - - - - - - - -
    **
    **  Quick observed place to CIRS, given the star-independent astrometry
    **  parameters.
    **
    **  Use of this function is appropriate when efficiency is important and
    **  where many star positions are all to be transformed for one date.
    **  The star-independent astrometry parameters can be obtained by
    **  calling iauApio[13] or iauApco[13].
    **
    **  Status:  support function.
    **
    **  Given:
    **     type   char[]     type of coordinates: "R", "H" or "A" (Note 1)
    **     ob1    double     observed Az, HA or RA (radians; Az is N=0,E=90)
    **     ob2    double     observed ZD or Dec (radians)
    **     astrom iauASTROM* star-independent astrometry parameters:
    **      pmt    double       PM time interval (SSB, Julian years)
    **      eb     double[3]    SSB to observer (vector, au)
    **      eh     double[3]    Sun to observer (unit vector)
    **      em     double       distance from Sun to observer (au)
    **      v      double[3]    barycentric observer velocity (vector, c)
    **      bm1    double       sqrt(1-|v|^2): reciprocal of Lorenz factor
    **      bpn    double[3][3] bias-precession-nutation matrix
    **      along  double       longitude + s' (radians)
    **      xpl    double       polar motion xp wrt local meridian (radians)
    **      ypl    double       polar motion yp wrt local meridian (radians)
    **      sphi   double       sine of geodetic latitude
    **      cphi   double       cosine of geodetic latitude
    **      diurab double       magnitude of diurnal aberration vector
    **      eral   double       "local" Earth rotation angle (radians)
    **      refa   double       refraction constant A (radians)
    **      refb   double       refraction constant B (radians)
    **
    **  Returned:
    **     ri     double*    CIRS right ascension (CIO-based, radians)
    **     di     double*    CIRS declination (radians)
    **
    **  Notes:
    **
    **  1) "Observed" Az,El means the position that would be seen by a
    **     perfect geodetically aligned theodolite.  This is related to
    **     the observed HA,Dec via the standard rotation, using the geodetic
    **     latitude (corrected for polar motion), while the observed HA and
    **     RA are related simply through the Earth rotation angle and the
    **     site longitude.  "Observed" RA,Dec or HA,Dec thus means the
    **     position that would be seen by a perfect equatorial with its
    **     polar axis aligned to the Earth's axis of rotation.  By removing
    **     from the observed place the effects of atmospheric refraction and
    **     diurnal aberration, the CIRS RA,Dec is obtained.
    **
    **  2) Only the first character of the type argument is significant.
    **     "R" or "r" indicates that ob1 and ob2 are the observed right
    **     ascension and declination;  "H" or "h" indicates that they are
    **     hour angle (west +ve) and declination;  anything else ("A" or
    **     "a" is recommended) indicates that ob1 and ob2 are azimuth (north
    **     zero, east 90 deg) and zenith distance.  (Zenith distance is used
    **     rather than altitude in order to reflect the fact that no
    **     allowance is made for depression of the horizon.)
    **
    **  3) The accuracy of the result is limited by the corrections for
    **     refraction, which use a simple A*tan(z) + B*tan^3(z) model.
    **     Providing the meteorological parameters are known accurately and
    **     there are no gross local effects, the predicted observed
    **     coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    **     (radio) for a zenith distance of less than 70 degrees, better
    **     than 30 arcsec (optical or radio) at 85 degrees and better than
    **     20 arcmin (optical) or 30 arcmin (radio) at the horizon.
    **
    **     Without refraction, the complementary functions iauAtioq and
    **     iauAtoiq are self-consistent to better than 1 microarcsecond all
    **     over the celestial sphere.  With refraction included, consistency
    **     falls off at high zenith distances, but is still better than
    **     0.05 arcsec at 85 degrees.
    **
    **  4) It is advisable to take great care with units, as even unlikely
    **     values of the input parameters are accepted and processed in
    **     accordance with the models used.
    **
    **  Called:
    **     iauS2c       spherical coordinates to unit vector
    **     iauC2s       p-vector to spherical
    **     iauAnp       normalize angle into range 0 to 2pi
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ri = c_double()
    di = c_double()
    strr = strtype.encode()
    _sofa.iauAtoiq(strr,ob1,ob2,astrom, byref(ri),byref(di))
    return ri.value, di.value

# iauLd
try:
    _sofa.iauLd.argtypes = [c_double, #bm 
                            ndpointer(shape=3, dtype=float), #p[3]
                            ndpointer(shape=3, dtype=float), #q[3]
                            ndpointer(shape=3, dtype=float), #e[3]
                            c_double, #em
                            c_double, #dlim
                            ndpointer(shape=3, dtype=float)]#p1
except AttributeError:
    pass  
def ld(bm,p,q,e,em,dlim):
    """
    

    Parameters
    ----------
    bm : double
        mass of the gravitating body (solar masses)
    p : double[3] 
        direction from observer to source (unit vector)
    q : double[3] 
        direction from body to source (unit vector).
    e : double[3] 
        direction from body to observer (unit vector).
    em : double
        distance from body to observer (au).
    dlim : double
        deflection limiter (Note 4).

    Returns
    -------
    p1 : double[3] 
        observer to deflected source (unit vector)
    /*
    **  - - - - - -
    **   i a u L d
    **  - - - - - -
    **
    **  Apply light deflection by a solar-system body, as part of
    **  transforming coordinate direction into natural direction.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     bm     double     mass of the gravitating body (solar masses)
    **     p      double[3]  direction from observer to source (unit vector)
    **     q      double[3]  direction from body to source (unit vector)
    **     e      double[3]  direction from body to observer (unit vector)
    **     em     double     distance from body to observer (au)
    **     dlim   double     deflection limiter (Note 4)
    **
    **  Returned:
    **     p1     double[3]  observer to deflected source (unit vector)
    **
    **  Notes:
    **
    **  1) The algorithm is based on Expr. (70) in Klioner (2003) and
    **     Expr. (7.63) in the Explanatory Supplement (Urban & Seidelmann
    **     2013), with some rearrangement to minimize the effects of machine
    **     precision.
    **
    **  2) The mass parameter bm can, as required, be adjusted in order to
    **     allow for such effects as quadrupole field.
    **
    **  3) The barycentric position of the deflecting body should ideally
    **     correspond to the time of closest approach of the light ray to
    **     the body.
    **
    **  4) The deflection limiter parameter dlim is phi^2/2, where phi is
    **     the angular separation (in radians) between source and body at
    **     which limiting is applied.  As phi shrinks below the chosen
    **     threshold, the deflection is artificially reduced, reaching zero
    **     for phi = 0.
    **
    **  5) The returned vector p1 is not normalized, but the consequential
    **     departure from unit magnitude is always negligible.
    **
    **  6) The arguments p and p1 can be the same array.
    **
    **  7) To accumulate total light deflection taking into account the
    **     contributions from several bodies, call the present function for
    **     each body in succession, in decreasing order of distance from the
    **     observer.
    **
    **  8) For efficiency, validation is omitted.  The supplied vectors must
    **     be of unit magnitude, and the deflection limiter non-zero and
    **     positive.
    **
    **  References:
    **
    **     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    **     the Astronomical Almanac, 3rd ed., University Science Books
    **     (2013).
    **
    **     Klioner, Sergei A., "A practical relativistic model for micro-
    **     arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).
    **
    **  Called:
    **     iauPdp       scalar product of two p-vectors
    **     iauPxp       vector product of two p-vectors
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    p1 = zeros(shape=3, dtype=float)
    _sofa.iauLd(bm,p,q, e,em,dlim,p1)
    return p1

#iauLdsun
try:
    _sofa.iauLdsun.argtypes = [ndpointer(shape=3, dtype=float), #p[3]
                               ndpointer(shape=3, dtype=float), #e[3]
                               c_double, #em
                               ndpointer(shape=3, dtype=float)]#p1
except AttributeError:
    pass  
def ldsun(p,e,em):
    """
    

    Parameters
    ----------
    p : double[3]
        direction from observer to star (unit vector)
    e : double[3]
        direction from Sun to observer (unit vector)
    em : double[3]
        distance from Sun to observer (au).

    Returns
    -------
    p1 : double[3]
        observer to deflected star (unit vector).

    **  Notes:
    **
    **  1) The source is presumed to be sufficiently distant that its
    **     directions seen from the Sun and the observer are essentially
    **     the same.
    **
    **  2) The deflection is restrained when the angle between the star and
    **     the center of the Sun is less than a threshold value, falling to
    **     zero deflection for zero separation.  The chosen threshold value
    **     is within the solar limb for all solar-system applications, and
    **     is about 5 arcminutes for the case of a terrestrial observer.
    **
    **  3) The arguments p and p1 can be the same array.
    **
    **  Called:
    **     iauLd        light deflection by a solar-system body
    **
    **  This revision:   2016 June 16
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    
    /*---------------------------
    """
    p1 = zeros(shape=3, dtype=float)
    _sofa.iauLdsun(p, e,em,p1)
    return p1
    

#iauPmpx
try:
    _sofa.iauPmpx.argtypes = [
        c_double, #rc
        c_double, #dc
        c_double, #pr
        c_double, #pd
        c_double, #px
        c_double, #rv,
        c_double, #pmt
        ndpointer(shape=3, dtype=float), #pob[3]
        ndpointer(shape=3, dtype=float)]#pco[3]]
except AttributeError:
    pass  
def pmpx(rc,dc,pr,pd,px,rv,pmt,pob):
    """
    

    Parameters
    ----------
    rc : double
        ICRS RA,Dec at catalog epoch (radians).
    dc : double
        ICRS RA,Dec at catalog epoch (radians).
    pr : double
        RA proper motion (radians/year; Note 1).
    pd : double
        Dec proper motion (radians/year)
    px : double
         parallax (arcsec).
    rv : double
        radial velocity (km/s, +ve if receding).
    pmt : double
        proper motion time interval (SSB, Julian years).
    pob : double[3] 
        DESCRIPTION.

    Returns
    -------
    pco : double[3]  
        coordinate direction (BCRS unit vector).
    /*
    **  - - - - - - - -
    **   i a u P m p x
    **  - - - - - - - -
    **
    **  Proper motion and parallax.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     rc,dc  double     ICRS RA,Dec at catalog epoch (radians)
    **     pr     double     RA proper motion (radians/year; Note 1)
    **     pd     double     Dec proper motion (radians/year)
    **     px     double     parallax (arcsec)
    **     rv     double     radial velocity (km/s, +ve if receding)
    **     pmt    double     proper motion time interval (SSB, Julian years)
    **     pob    double[3]  SSB to observer vector (au)
    **
    **  Returned:
    **     pco    double[3]  coordinate direction (BCRS unit vector)
    **
    **  Notes:
    **
    **  1) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
    **
    **  2) The proper motion time interval is for when the starlight
    **     reaches the solar system barycenter.
    **
    **  3) To avoid the need for iteration, the Roemer effect (i.e. the
    **     small annual modulation of the proper motion coming from the
    **     changing light time) is applied approximately, using the
    **     direction of the star at the catalog epoch.
    **
    **  References:
    **
    **     1984 Astronomical Almanac, pp B39-B41.
    **
    **     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    **     the Astronomical Almanac, 3rd ed., University Science Books
    **     (2013), Section 7.2.
    **
    **  Called:
    **     iauPdp       scalar product of two p-vectors
    **     iauPn        decompose p-vector into modulus and direction
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    pco = zeros(shape=3, dtype=float)
    _sofa.iauPmpx(rc,dc,pr,pd,px,rv,pmt,pob,pco)
    return pco


#iauPmsafe
pmsafe_msg = {0: 'OK', # Unused
                        -1 : 'system error (should not occur)',
                        1 : 'distance overridden (Note 6)',
                        2 : 'excessive velocity (Note 7)',
                        4 : 'solution didnot converge (Note 8)'}
try:
    _sofa.iauPmsafe.argtypes = [
    c_double,  # ra1 \RA(radians), before
    c_double,  # dec1 \DEC(radians), before
    c_double,  # pmr1 \RA proper motion(radians/year), before
    c_double,  # pmd1 \DEC proper motion(radians/year), before
    c_double,  # px1 \parallax(arcsec), before
    c_double,  # rv1 \radial velocity(km/s,positive if receding), before
    c_double,  # ep1a \"before" epoch, part A
    c_double,  # ep1b \"before" epoch, part B
    c_double,  # ep2a \"after" epoch, part A
    c_double,  # ep2b \"after" epoch, part B
    POINTER(c_double),  # ra2 \ right ascension (radians), after
    POINTER(c_double),  # dec2 \ declination (radians), after
    POINTER(c_double),  # pmr2 \ RA proper motion (radians/year), after
    POINTER(c_double),  # pmd2 \ Dec proper motion (radians/year), after
    POINTER(c_double),  # px2 \ parallax (arcseconds), after
    POINTER(c_double),  # rv2 \ radial velocity (km/s, +ve = receding), after
]
except AttributeError:
    pass  
def pmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b):
    """
    

    Parameters
    ----------
    ra1 : TYPE
        DESCRIPTION.
    dec1 : TYPE
        DESCRIPTION.
    pmr1 : TYPE
        DESCRIPTION.
    pmd1 : TYPE
        DESCRIPTION.
    px1 : TYPE
        DESCRIPTION.
    rv1 : TYPE
        DESCRIPTION.
    ep1a : TYPE
        DESCRIPTION.
    ep1b : TYPE
        DESCRIPTION.
    ep2a : TYPE
        DESCRIPTION.
    ep2b : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.
    /*
    **  - - - - - - - - - -
    **   i a u P m s a f e
    **  - - - - - - - - - -
    **
    **  Star proper motion:  update star catalog data for space motion, with
    **  special handling to handle the zero parallax case.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ra1    double      right ascension (radians), before
    **     dec1   double      declination (radians), before
    **     pmr1   double      RA proper motion (radians/year), before
    **     pmd1   double      Dec proper motion (radians/year), before
    **     px1    double      parallax (arcseconds), before
    **     rv1    double      radial velocity (km/s, +ve = receding), before
    **     ep1a   double      "before" epoch, part A (Note 1)
    **     ep1b   double      "before" epoch, part B (Note 1)
    **     ep2a   double      "after" epoch, part A (Note 1)
    **     ep2b   double      "after" epoch, part B (Note 1)
    **
    **  Returned:
    **     ra2    double      right ascension (radians), after
    **     dec2   double      declination (radians), after
    **     pmr2   double      RA proper motion (radians/year), after
    **     pmd2   double      Dec proper motion (radians/year), after
    **     px2    double      parallax (arcseconds), after
    **     rv2    double      radial velocity (km/s, +ve = receding), after
    **
    **  Returned (function value):
    **            int         status:
    **                         -1 = system error (should not occur)
    **                          0 = no warnings or errors
    **                          1 = distance overridden (Note 6)
    **                          2 = excessive velocity (Note 7)
    **                          4 = solution didn't converge (Note 8)
    **                       else = binary logical OR of the above warnings
    **
    **  Notes:
    **
    **  1) The starting and ending TDB epochs ep1a+ep1b and ep2a+ep2b are
    **     Julian Dates, apportioned in any convenient way between the two
    **     parts (A and B).  For example, JD(TDB)=2450123.7 could be
    **     expressed in any of these ways, among others:
    **
    **            epNa            epNb
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in cases
    **     where the loss of several decimal digits of resolution is
    **     acceptable.  The J2000 method is best matched to the way the
    **     argument is handled internally and will deliver the optimum
    **     resolution.  The MJD method and the date & time methods are both
    **     good compromises between resolution and convenience.
    **
    **  2) In accordance with normal star-catalog conventions, the object's
    **     right ascension and declination are freed from the effects of
    **     secular aberration.  The frame, which is aligned to the catalog
    **     equator and equinox, is Lorentzian and centered on the SSB.
    **
    **     The proper motions are the rate of change of the right ascension
    **     and declination at the catalog epoch and are in radians per TDB
    **     Julian year.
    **
    **     The parallax and radial velocity are in the same frame.
    **
    **  3) Care is needed with units.  The star coordinates are in radians
    **     and the proper motions in radians per Julian year, but the
    **     parallax is in arcseconds.
    **
    **  4) The RA proper motion is in terms of coordinate angle, not true
    **     angle.  If the catalog uses arcseconds for both RA and Dec proper
    **     motions, the RA proper motion will need to be divided by cos(Dec)
    **     before use.
    **
    **  5) Straight-line motion at constant speed, in the inertial frame, is
    **     assumed.
    **
    **  6) An extremely small (or zero or negative) parallax is overridden
    **     to ensure that the object is at a finite but very large distance,
    **     but not so large that the proper motion is equivalent to a large
    **     but safe speed (about 0.1c using the chosen constant).  A warning
    **     status of 1 is added to the status if this action has been taken.
    **
    **  7) If the space velocity is a significant fraction of c (see the
    **     constant VMAX in the function iauStarpv), it is arbitrarily set
    **     to zero.  When this action occurs, 2 is added to the status.
    **
    **  8) The relativistic adjustment carried out in the iauStarpv function
    **     involves an iterative calculation.  If the process fails to
    **     converge within a set number of iterations, 4 is added to the
    **     status.
    **
    **  Called:
    **     iauSeps      angle between two points
    **     iauStarpm    update star catalog data for space motion
    **
    **  This revision:   2014 July 1
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ra2 = c_double()
    dec2 = c_double()
    pmr2 = c_double()
    pmd2 = c_double()
    px2 = c_double()
    rv2 = c_double()
    s = _sofa.iauPmsafe(ra1, dec1, pmr1, pmd1, px1, rv1, ep1a, ep1b, ep2a, ep2b,
                    byref(ra2), byref(dec2), byref(pmr2), byref(pmd2), byref(px2), byref(rv2))
    if s != 0:
        warnings.warn(pmsafe_msg[s], UserWarning, 2)
    return ra2.value, dec2.value, pmr2.value, pmd2.value, px2.value, rv2.value


#iauPvtob
try:
    _sofa.iauPvtob.argtypes = [
    c_double,  # elong 
    c_double,  # phi
    c_double,  # hm
    c_double,  # xp
    c_double,  # yp
    c_double,  # sp
    c_double,  # theta
    ndpointer(shape = (2,3),dtype = float),  #pv[2][3]
]
except AttributeError:
    pass  
def pvtob(elong,phi,hm,xp,yp,sp,theta):
    """
    

    Parameters
    ----------
    elong : double
        longitude (radians, east +ve, Note 1).
    phi : TYPE
        latitude (geodetic, radians, Note 1).
    hm : TYPE
        height above ref. ellipsoid (geodetic, m).
    xp : TYPE
        coordinates of the pole (radians, Note 2).
    yp : TYPE
        DESCRIPTION.
    sp : TYPE
        DESCRIPTION.
    theta : TYPE
        the TIO locator s' (radians, Note 2).

    Returns
    -------
    pv[2][3] : TYPE
        double[2][3] position/velocity vector (m, m/s, CIRS).
        /*
    **  - - - - - - - - -
    **   i a u P v t o b
    **  - - - - - - - - -
    **
    **  Position and velocity of a terrestrial observing station.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     elong   double       longitude (radians, east +ve, Note 1)
    **     phi     double       latitude (geodetic, radians, Note 1)
    **     hm      double       height above ref. ellipsoid (geodetic, m)
    **     xp,yp   double       coordinates of the pole (radians, Note 2)
    **     sp      double       the TIO locator s' (radians, Note 2)
    **     theta   double       Earth rotation angle (radians, Note 3)
    **
    **  Returned:
    **     pv      double[2][3] position/velocity vector (m, m/s, CIRS)
    **
    **  Notes:
    **
    **  1) The terrestrial coordinates are with respect to the WGS84
    **     reference ellipsoid.
    **
    **  2) xp and yp are the coordinates (in radians) of the Celestial
    **     Intermediate Pole with respect to the International Terrestrial
    **     Reference System (see IERS Conventions), measured along the
    **     meridians 0 and 90 deg west respectively.  sp is the TIO locator
    **     s', in radians, which positions the Terrestrial Intermediate
    **     Origin on the equator.  For many applications, xp, yp and
    **     (especially) sp can be set to zero.
    **
    **  3) If theta is Greenwich apparent sidereal time instead of Earth
    **     rotation angle, the result is with respect to the true equator
    **     and equinox of date, i.e. with the x-axis at the equinox rather
    **     than the celestial intermediate origin.
    **
    **  4) The velocity units are meters per UT1 second, not per SI second.
    **     This is unlikely to have any practical consequences in the modern
    **     era.
    **
    **  5) No validation is performed on the arguments.  Error cases that
    **     could lead to arithmetic exceptions are trapped by the iauGd2gc
    **     function, and the result set to zeros.
    **
    **  References:
    **
    **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    **     IERS Technical Note No. 32, BKG (2004)
    **
    **     Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    **     the Astronomical Almanac, 3rd ed., University Science Books
    **     (2013), Section 7.4.3.3.
    **
    **  Called:
    **     iauGd2gc     geodetic to geocentric transformation
    **     iauPom00     polar motion matrix
    **     iauTrxp      product of transpose of r-matrix and p-vector
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    pv = zeros(shape = (2,3), dtype = float)
    _sofa.iauPvtob(elong,phi,hm,xp,yp,sp,theta,pv)
    return pv


#iauPvstar
pvstar_msg = {0: 'OK', # Unused
                -1:'superluminal speed (Note 5)',
                -2:'null position vector'
                }
try:
    _sofa.iauPvstar.argtypes = [
    ndpointer(shape = (2,3), dtype = float),  # pv 
    POINTER(c_double),  # ra
    POINTER(c_double),  # dec
    POINTER(c_double),  # pmr
    POINTER(c_double),  # pmd
    POINTER(c_double),  # px
    POINTER(c_double)   # rv
]
except AttributeError:
    pass  
def pvstar(pv):
    """
    

    Parameters
    ----------
    pv : double[2][3]
         pv-vector (au, au/day).

    Returns
    -------
    ra : double
        right ascension (radians)
    dec : double
        declination (radians)
    pmr : double
        RA proper motion (radians/year).
    pmd : double
        Dec proper motion (radians/year).
    px : double
        parallax (arcsec).
    rv : double
        radial velocity (km/s, positive = receding).
    /*
    **  - - - - - - - - - -
    **   i a u P v s t a r
    **  - - - - - - - - - -
    **
    **  Convert star position+velocity vector to catalog coordinates.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given (Note 1):
    **     pv     double[2][3]   pv-vector (au, au/day)
    **
    **  Returned (Note 2):
    **     ra     double         right ascension (radians)
    **     dec    double         declination (radians)
    **     pmr    double         RA proper motion (radians/year)
    **     pmd    double         Dec proper motion (radians/year)
    **     px     double         parallax (arcsec)
    **     rv     double         radial velocity (km/s, positive = receding)
    **
    **  Returned (function value):
    **            int            status:
    **                              0 = OK
    **                             -1 = superluminal speed (Note 5)
    **                             -2 = null position vector
    **
    **  Notes:
    **
    **  1) The specified pv-vector is the coordinate direction (and its rate
    **     of change) for the date at which the light leaving the star
    **     reached the solar-system barycenter.
    **
    **  2) The star data returned by this function are "observables" for an
    **     imaginary observer at the solar-system barycenter.  Proper motion
    **     and radial velocity are, strictly, in terms of barycentric
    **     coordinate time, TCB.  For most practical applications, it is
    **     permissible to neglect the distinction between TCB and ordinary
    **     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
    **     limited by the intrinsic accuracy of the proper-motion and
    **     radial-velocity data;  moreover, the supplied pv-vector is likely
    **     to be merely an intermediate result (for example generated by the
    **     function iauStarpv), so that a change of time unit will cancel
    **     out overall.
    **
    **     In accordance with normal star-catalog conventions, the object's
    **     right ascension and declination are freed from the effects of
    **     secular aberration.  The frame, which is aligned to the catalog
    **     equator and equinox, is Lorentzian and centered on the SSB.
    **
    **     Summarizing, the specified pv-vector is for most stars almost
    **     identical to the result of applying the standard geometrical
    **     "space motion" transformation to the catalog data.  The
    **     differences, which are the subject of the Stumpff paper cited
    **     below, are:
    **
    **     (i) In stars with significant radial velocity and proper motion,
    **     the constantly changing light-time distorts the apparent proper
    **     motion.  Note that this is a classical, not a relativistic,
    **     effect.
    **
    **     (ii) The transformation complies with special relativity.
    **
    **  3) Care is needed with units.  The star coordinates are in radians
    **     and the proper motions in radians per Julian year, but the
    **     parallax is in arcseconds; the radial velocity is in km/s, but
    **     the pv-vector result is in au and au/day.
    **
    **  4) The proper motions are the rate of change of the right ascension
    **     and declination at the catalog epoch and are in radians per Julian
    **     year.  The RA proper motion is in terms of coordinate angle, not
    **     true angle, and will thus be numerically larger at high
    **     declinations.
    **
    **  5) Straight-line motion at constant speed in the inertial frame is
    **     assumed.  If the speed is greater than or equal to the speed of
    **     light, the function aborts with an error status.
    **
    **  6) The inverse transformation is performed by the function iauStarpv.
    **
    **  Called:
    **     iauPn        decompose p-vector into modulus and direction
    **     iauPdp       scalar product of two p-vectors
    **     iauSxp       multiply p-vector by scalar
    **     iauPmp       p-vector minus p-vector
    **     iauPm        modulus of p-vector
    **     iauPpp       p-vector plus p-vector
    **     iauPv2s      pv-vector to spherical
    **     iauAnp       normalize angle into range 0 to 2pi
    **
    **  Reference:
    **
    **     Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
    **
    **  This revision:  2017 March 16
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ra = c_double()
    dec = c_double()
    pmr = c_double()
    pmd = c_double()
    px = c_double()
    rv = c_double()
    s = _sofa.iauPvstar(pv,byref(ra),byref(dec),byref(pmr),byref(pmd),byref(px),byref(rv))
    if s != 0:
        warnings.warn(pvstar_msg[s], UserWarning, 2)
    return ra.value, dec.value, pmr.value, pmd.value, px.value, rv.value


# iauRefco
try:
    _sofa.iauRefco.argtypes = [
    c_double,  # phpa
    c_double,  # tc
    c_double,  # rh
    c_double,  # wl
    POINTER(c_double),  # refa
    POINTER(c_double)   # refb
]
except AttributeError:
    pass  
def refco(phpa, tc, rh, wl):
    """
    

    Parameters
    ----------
    phpa : double
        pressure at the observer (hPa = millibar).
    tc : double
        ambient temperature at the observer (deg C).
    rh : double
        relative humidity at the observer (range 0-1).
    wl : double
        wavelength (micrometers).

    Returns
    -------
    refa : double
        tan Z coefficient (radians).
    refb : double
        tan^3 Z coefficient (radians).
    /*
    **  - - - - - - - - -
    **   i a u R e f c o
    **  - - - - - - - - -
    **
    **  Determine the constants A and B in the atmospheric refraction model
    **  dZ = A tan Z + B tan^3 Z.
    **
    **  Z is the "observed" zenith distance (i.e. affected by refraction)
    **  and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
    **  zenith distance.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **    phpa   double    pressure at the observer (hPa = millibar)
    **    tc     double    ambient temperature at the observer (deg C)
    **    rh     double    relative humidity at the observer (range 0-1)
    **    wl     double    wavelength (micrometers)
    **
    **  Returned:
    **    refa   double*   tan Z coefficient (radians)
    **    refb   double*   tan^3 Z coefficient (radians)
    **
    **  Notes:
    **
    **  1) The model balances speed and accuracy to give good results in
    **     applications where performance at low altitudes is not paramount.
    **     Performance is maintained across a range of conditions, and
    **     applies to both optical/IR and radio.
    **
    **  2) The model omits the effects of (i) height above sea level (apart
    **     from the reduced pressure itself), (ii) latitude (i.e. the
    **     flattening of the Earth), (iii) variations in tropospheric lapse
    **     rate and (iv) dispersive effects in the radio.
    **
    **     The model was tested using the following range of conditions:
    **
    **       lapse rates 0.0055, 0.0065, 0.0075 deg/meter
    **       latitudes 0, 25, 50, 75 degrees
    **       heights 0, 2500, 5000 meters ASL
    **       pressures mean for height -10% to +5% in steps of 5%
    **       temperatures -10 deg to +20 deg with respect to 280 deg at SL
    **       relative humidity 0, 0.5, 1
    **       wavelengths 0.4, 0.6, ... 2 micron, + radio
    **       zenith distances 15, 45, 75 degrees
    **
    **     The accuracy with respect to raytracing through a model
    **     atmosphere was as follows:
    **
    **                            worst         RMS
    **
    **       optical/IR           62 mas       8 mas
    **       radio               319 mas      49 mas
    **
    **     For this particular set of conditions:
    **
    **       lapse rate 0.0065 K/meter
    **       latitude 50 degrees
    **       sea level
    **       pressure 1005 mb
    **       temperature 280.15 K
    **       humidity 80%
    **       wavelength 5740 Angstroms
    **
    **     the results were as follows:
    **
    **       ZD       raytrace     iauRefco   Saastamoinen
    **
    **       10         10.27        10.27        10.27
    **       20         21.19        21.20        21.19
    **       30         33.61        33.61        33.60
    **       40         48.82        48.83        48.81
    **       45         58.16        58.18        58.16
    **       50         69.28        69.30        69.27
    **       55         82.97        82.99        82.95
    **       60        100.51       100.54       100.50
    **       65        124.23       124.26       124.20
    **       70        158.63       158.68       158.61
    **       72        177.32       177.37       177.31
    **       74        200.35       200.38       200.32
    **       76        229.45       229.43       229.42
    **       78        267.44       267.29       267.41
    **       80        319.13       318.55       319.10
    **
    **      deg        arcsec       arcsec       arcsec
    **
    **     The values for Saastamoinen's formula (which includes terms
    **     up to tan^5) are taken from Hohenkerk and Sinclair (1985).
    **
    **  3) A wl value in the range 0-100 selects the optical/IR case and is
    **     wavelength in micrometers.  Any value outside this range selects
    **     the radio case.
    **
    **  4) Outlandish input parameters are silently limited to
    **     mathematically safe values.  Zero pressure is permissible, and
    **     causes zeroes to be returned.
    **
    **  5) The algorithm draws on several sources, as follows:
    **
    **     a) The formula for the saturation vapour pressure of water as
    **        a function of temperature and temperature is taken from
    **        Equations (A4.5-A4.7) of Gill (1982).
    **
    **     b) The formula for the water vapour pressure, given the
    **        saturation pressure and the relative humidity, is from
    **        Crane (1976), Equation (2.5.5).
    **
    **     c) The refractivity of air is a function of temperature,
    **        total pressure, water-vapour pressure and, in the case
    **        of optical/IR, wavelength.  The formulae for the two cases are
    **        developed from Hohenkerk & Sinclair (1985) and Rueger (2002).
    **
    **     d) The formula for beta, the ratio of the scale height of the
    **        atmosphere to the geocentric distance of the observer, is
    **        an adaption of Equation (9) from Stone (1996).  The
    **        adaptations, arrived at empirically, consist of (i) a small
    **        adjustment to the coefficient and (ii) a humidity term for the
    **        radio case only.
    **
    **     e) The formulae for the refraction constants as a function of
    **        n-1 and beta are from Green (1987), Equation (4.31).
    **
    **  References:
    **
    **     Crane, R.K., Meeks, M.L. (ed), "Refraction Effects in the Neutral
    **     Atmosphere", Methods of Experimental Physics: Astrophysics 12B,
    **     Academic Press, 1976.
    **
    **     Gill, Adrian E., "Atmosphere-Ocean Dynamics", Academic Press,
    **     1982.
    **
    **     Green, R.M., "Spherical Astronomy", Cambridge University Press,
    **     1987.
    **
    **     Hohenkerk, C.Y., & Sinclair, A.T., NAO Technical Note No. 63,
    **     1985.
    **
    **     Rueger, J.M., "Refractive Index Formulae for Electronic Distance
    **     Measurement with Radio and Millimetre Waves", in Unisurv Report
    **     S-68, School of Surveying and Spatial Information Systems,
    **     University of New South Wales, Sydney, Australia, 2002.
    **
    **     Stone, Ronald C., P.A.S.P. 108, 1051-1058, 1996.
    **
    **  This revision:   2013 October 9
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    refa = c_double()
    refb = c_double()
    _sofa.iauRefco(phpa, tc, rh,wl, byref(refa), byref(refb))
    return refa.value, refb.value


# iauStarpm
starpm_msg = {0: 'no warnings or errors', 
                -1: 'ssystem error (should not occur)',
                1: 'distance overridden (Note 6)',
                2: 'excessive velocity (Note 7)',
                4: 'solution didnot converge (Note 8)'
                }
try:
    _sofa.iauStarpm.argtypes = [
    c_double,  # ra1
    c_double,  # dec1
    c_double,  # pmr1
    c_double,  # pmd1
    c_double,  # px1
    c_double,  # rv1
    c_double,  #ep1a
    c_double,  #ep1b
    c_double,  #ep2a
    c_double,  #ep2b
    POINTER(c_double),  # ra2
    POINTER(c_double),  # dec2
    POINTER(c_double),  # pmr2
    POINTER(c_double),  # pmd2
    POINTER(c_double),  # px2
    POINTER(c_double),  # rv2
]
except AttributeError:
    pass  
def starpm(ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b):
    """
    Star proper motion:  update star catalog data for space motion.

    Parameters
    ----------
    ra1 : TYPE
        right ascension (radians), before.
    dec1 : TYPE
        declination (radians), before.
    pmr1 : TYPE
        RA proper motion (radians/year), before.
    pmd1 : TYPE
        Dec proper motion (radians/year), before.
    px1 : TYPE
        parallax (arcseconds), before.
    rv1 : TYPE
        radial velocity (km/s, +ve = receding), before.
    ep1a : TYPE
        "before" epoch, part A (Note 1).
    ep1b : TYPE
        "before" epoch, part B (Note 1).
    ep2a : TYPE
        "after" epoch, part A (Note 1).
    ep2b : TYPE
        "after" epoch, part B (Note 1).

    Returns
    -------
    ra2 : TYPE
        right ascension (radians), after.
    dec2 : 
        declination (radians), after.
    pmr2 :
        RA proper motion (radians/year), after.
    pmd2 :
        Dec proper motion (radians/year), after.
    px2 : 
        parallax (arcseconds), after.
    rv2 : 
        radial velocity (km/s, +ve = receding), after.
    /*
    **  - - - - - - - - - -
    **   i a u S t a r p m
    **  - - - - - - - - - -
    **
    **  Star proper motion:  update star catalog data for space motion.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     ra1    double     right ascension (radians), before
    **     dec1   double     declination (radians), before
    **     pmr1   double     RA proper motion (radians/year), before
    **     pmd1   double     Dec proper motion (radians/year), before
    **     px1    double     parallax (arcseconds), before
    **     rv1    double     radial velocity (km/s, +ve = receding), before
    **     ep1a   double     "before" epoch, part A (Note 1)
    **     ep1b   double     "before" epoch, part B (Note 1)
    **     ep2a   double     "after" epoch, part A (Note 1)
    **     ep2b   double     "after" epoch, part B (Note 1)
    **
    **  Returned:
    **     ra2    double     right ascension (radians), after
    **     dec2   double     declination (radians), after
    **     pmr2   double     RA proper motion (radians/year), after
    **     pmd2   double     Dec proper motion (radians/year), after
    **     px2    double     parallax (arcseconds), after
    **     rv2    double     radial velocity (km/s, +ve = receding), after
    **
    **  Returned (function value):
    **            int        status:
    **                          -1 = system error (should not occur)
    **                           0 = no warnings or errors
    **                           1 = distance overridden (Note 6)
    **                           2 = excessive velocity (Note 7)
    **                           4 = solution didn't converge (Note 8)
    **                        else = binary logical OR of the above warnings
    **
    **  Notes:
    **
    **  1) The starting and ending TDB dates ep1a+ep1b and ep2a+ep2b are
    **     Julian Dates, apportioned in any convenient way between the two
    **     parts (A and B).  For example, JD(TDB)=2450123.7 could be
    **     expressed in any of these ways, among others:
    **
    **             epna          epnb
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable.  The J2000 method is best matched to the way
    **     the argument is handled internally and will deliver the
    **     optimum resolution.  The MJD method and the date & time methods
    **     are both good compromises between resolution and convenience.
    **
    **  2) In accordance with normal star-catalog conventions, the object's
    **     right ascension and declination are freed from the effects of
    **     secular aberration.  The frame, which is aligned to the catalog
    **     equator and equinox, is Lorentzian and centered on the SSB.
    **
    **     The proper motions are the rate of change of the right ascension
    **     and declination at the catalog epoch and are in radians per TDB
    **     Julian year.
    **
    **     The parallax and radial velocity are in the same frame.
    **
    **  3) Care is needed with units.  The star coordinates are in radians
    **     and the proper motions in radians per Julian year, but the
    **     parallax is in arcseconds.
    **
    **  4) The RA proper motion is in terms of coordinate angle, not true
    **     angle.  If the catalog uses arcseconds for both RA and Dec proper
    **     motions, the RA proper motion will need to be divided by cos(Dec)
    **     before use.
    **
    **  5) Straight-line motion at constant speed, in the inertial frame,
    **     is assumed.
    **
    **  6) An extremely small (or zero or negative) parallax is interpreted
    **     to mean that the object is on the "celestial sphere", the radius
    **     of which is an arbitrary (large) value (see the iauStarpv
    **     function for the value used).  When the distance is overridden in
    **     this way, the status, initially zero, has 1 added to it.
    **
    **  7) If the space velocity is a significant fraction of c (see the
    **     constant VMAX in the function iauStarpv), it is arbitrarily set
    **     to zero.  When this action occurs, 2 is added to the status.
    **
    **  8) The relativistic adjustment carried out in the iauStarpv function
    **     involves an iterative calculation.  If the process fails to
    **     converge within a set number of iterations, 4 is added to the
    **     status.
    **
    **  Called:
    **     iauStarpv    star catalog data to space motion pv-vector
    **     iauPvu       update a pv-vector
    **     iauPdp       scalar product of two p-vectors
    **     iauPvstar    space motion pv-vector to star catalog data
    **
    **  This revision:  2013 June 18
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    ra2 = c_double()
    dec2 = c_double()
    pmr2 = c_double()
    pmd2 = c_double()
    px2 = c_double()
    rv2 = c_double()
    s = _sofa.iauStarpm(ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b,byref(ra2),byref(dec2),byref(pmr2),byref(pmd2),byref(px2),byref(rv2))
    if s in starpm_msg:
        warnings.warn(starpm_msg[s], UserWarning, 2)
    else:
        warnings.warn('binary logical OR of the above warnings')
    return ra2.value, dec2.value, pmr2.value, pmd2.value, px2.value, rv2.value


# iauStarpv
starpv_msg = {0: 'no warnings or errors', 
                1: 'distance overridden (Note 6)',
                2: 'excessive speed (Note 7)',
                4: 'solution didnot converge (Note 8)'
                }
try:
    _sofa.iauStarpv.argtypes = [
    c_double,  # ra
    c_double,  # dec
    c_double,  # pmr
    c_double,  # pmd
    c_double,  # px
    c_double,  # rv
    ndpointer(shape = (2,3), dtype = float),  # pv[2][3]
]
except AttributeError:
    pass  
def starpv(ra,dec,pmr,pmd,px,rv):
    """
    Convert star catalog coordinates to position+velocity vector.

    Parameters
    ----------
    ra : double
        right ascension (radians).
    dec : double
        declination (radians).
    pmr : double
        RA proper motion (radians/year).
    pmd : double
        Dec proper motion (radians/year).
    px : double
        parallax (arcseconds).
    rv : double
        radial velocity (km/s, positive = receding).

    Returns
    -------
    pv : double[2][3]
        pv-vector (au, au/day).
    /*
    **  - - - - - - - - - -
    **   i a u S t a r p v
    **  - - - - - - - - - -
    **
    **  Convert star catalog coordinates to position+velocity vector.
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given (Note 1):
    **     ra     double        right ascension (radians)
    **     dec    double        declination (radians)
    **     pmr    double        RA proper motion (radians/year)
    **     pmd    double        Dec proper motion (radians/year)
    **     px     double        parallax (arcseconds)
    **     rv     double        radial velocity (km/s, positive = receding)
    **
    **  Returned (Note 2):
    **     pv     double[2][3]  pv-vector (au, au/day)
    **
    **  Returned (function value):
    **            int           status:
    **                              0 = no warnings
    **                              1 = distance overridden (Note 6)
    **                              2 = excessive speed (Note 7)
    **                              4 = solution didn't converge (Note 8)
    **                           else = binary logical OR of the above
    **
    **  Notes:
    **
    **  1) The star data accepted by this function are "observables" for an
    **     imaginary observer at the solar-system barycenter.  Proper motion
    **     and radial velocity are, strictly, in terms of barycentric
    **     coordinate time, TCB.  For most practical applications, it is
    **     permissible to neglect the distinction between TCB and ordinary
    **     "proper" time on Earth (TT/TAI).  The result will, as a rule, be
    **     limited by the intrinsic accuracy of the proper-motion and
    **     radial-velocity data;  moreover, the pv-vector is likely to be
    **     merely an intermediate result, so that a change of time unit
    **     would cancel out overall.
    **
    **     In accordance with normal star-catalog conventions, the object's
    **     right ascension and declination are freed from the effects of
    **     secular aberration.  The frame, which is aligned to the catalog
    **     equator and equinox, is Lorentzian and centered on the SSB.
    **
    **  2) The resulting position and velocity pv-vector is with respect to
    **     the same frame and, like the catalog coordinates, is freed from
    **     the effects of secular aberration.  Should the "coordinate
    **     direction", where the object was located at the catalog epoch, be
    **     required, it may be obtained by calculating the magnitude of the
    **     position vector pv[0][0-2] dividing by the speed of light in
    **     au/day to give the light-time, and then multiplying the space
    **     velocity pv[1][0-2] by this light-time and adding the result to
    **     pv[0][0-2].
    **
    **     Summarizing, the pv-vector returned is for most stars almost
    **     identical to the result of applying the standard geometrical
    **     "space motion" transformation.  The differences, which are the
    **     subject of the Stumpff paper referenced below, are:
    **
    **     (i) In stars with significant radial velocity and proper motion,
    **     the constantly changing light-time distorts the apparent proper
    **     motion.  Note that this is a classical, not a relativistic,
    **     effect.
    **
    **     (ii) The transformation complies with special relativity.
    **
    **  3) Care is needed with units.  The star coordinates are in radians
    **     and the proper motions in radians per Julian year, but the
    **     parallax is in arcseconds; the radial velocity is in km/s, but
    **     the pv-vector result is in au and au/day.
    **
    **  4) The RA proper motion is in terms of coordinate angle, not true
    **     angle.  If the catalog uses arcseconds for both RA and Dec proper
    **     motions, the RA proper motion will need to be divided by cos(Dec)
    **     before use.
    **
    **  5) Straight-line motion at constant speed, in the inertial frame,
    **     is assumed.
    **
    **  6) An extremely small (or zero or negative) parallax is interpreted
    **     to mean that the object is on the "celestial sphere", the radius
    **     of which is an arbitrary (large) value (see the constant PXMIN).
    **     When the distance is overridden in this way, the status,
    **     initially zero, has 1 added to it.
    **
    **  7) If the space velocity is a significant fraction of c (see the
    **     constant VMAX), it is arbitrarily set to zero.  When this action
    **     occurs, 2 is added to the status.
    **
    **  8) The relativistic adjustment involves an iterative calculation.
    **     If the process fails to converge within a set number (IMAX) of
    **     iterations, 4 is added to the status.
    **
    **  9) The inverse transformation is performed by the function
    **     iauPvstar.
    **
    **  Called:
    **     iauS2pv      spherical coordinates to pv-vector
    **     iauPm        modulus of p-vector
    **     iauZp        zero p-vector
    **     iauPn        decompose p-vector into modulus and direction
    **     iauPdp       scalar product of two p-vectors
    **     iauSxp       multiply p-vector by scalar
    **     iauPmp       p-vector minus p-vector
    **     iauPpp       p-vector plus p-vector
    **
    **  Reference:
    **
    **     Stumpff, P., 1985, Astron.Astrophys. 144, 232-240.
    **
    **  This revision:  2017 March 16
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    pv = zeros(shape = (2,3), dtype = float)
    s = _sofa.iauStarpv(ra,dec,pmr,pmd,px,rv,pv)
    if s in starpv_msg:
        warnings.warn(starpv_msg[s], UserWarning, 2)
    else:
        warnings.warn('binary logical OR of the above warnings')
    return pv 


##############################################################################
  #  Time scales
##############################################################################
# iauD2DTF
d2dtf_msg = {0: 'OK', # Unused
                -1:'unacceptable date (Note 6)',
                +1:'dubious year (Note 5)'}
_sofa.iauD2dtf.argtypes = [
    POINTER(c_char),   # scale, char[] 
    c_int, # ndp
    c_double, # d1
    c_double, # d2
    POINTER(c_int), #iy
    POINTER(c_int), #im
    POINTER(c_int), #id
    ndpointer(shape  =(4,1), dtype = int) #ihmsf[4]
]
def d2dtf(scale, ndp, d1, d2):
    """
    Format for output a 2-part Julian Date (or in the case of UTC a
 quasi-JD form that includes special provision for leap seconds).

    Parameters
    ----------
    scale : char[]
        time scale ID (Note 1).
    ndp : int
        resolution (Note 2).
    d1 : float
        time as a 2-part Julian Date (Notes 3,4).
    d2 : float
        time as a 2-part Julian Date (Notes 3,4).

    Returns
    -------
    iy : int
    im : int
    id :int
        year, month, day in Gregorian calendar (Note 5).
    ihmsf : int[4]
        hours, minutes, seconds, fraction (Note 1).
    /*
    **  - - - - - - - - -
    **   i a u D 2 d t f
    **  - - - - - - - - -
    **
    **  Format for output a 2-part Julian Date (or in the case of UTC a
    **  quasi-JD form that includes special provision for leap seconds).
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     scale     char[]  time scale ID (Note 1)
    **     ndp       int     resolution (Note 2)
    **     d1,d2     double  time as a 2-part Julian Date (Notes 3,4)
    **
    **  Returned:
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 5)
    **     ihmsf     int[4]  hours, minutes, seconds, fraction (Note 1)
    **
    **  Returned (function value):
    **               int     status: +1 = dubious year (Note 5)
    **                                0 = OK
    **                               -1 = unacceptable date (Note 6)
    **
    **  Notes:
    **
    **  1) scale identifies the time scale.  Only the value "UTC" (in upper
    **     case) is significant, and enables handling of leap seconds (see
    **     Note 4).
    **
    **  2) ndp is the number of decimal places in the seconds field, and can
    **     have negative as well as positive values, such as:
    **
    **     ndp         resolution
    **     -4            1 00 00
    **     -3            0 10 00
    **     -2            0 01 00
    **     -1            0 00 10
    **      0            0 00 01
    **      1            0 00 00.1
    **      2            0 00 00.01
    **      3            0 00 00.001
    **
    **     The limits are platform dependent, but a safe range is -5 to +9.
    **
    **  3) d1+d2 is Julian Date, apportioned in any convenient way between
    **     the two arguments, for example where d1 is the Julian Day Number
    **     and d2 is the fraction of a day.  In the case of UTC, where the
    **     use of JD is problematical, special conventions apply:  see the
    **     next note.
    **
    **  4) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The SOFA internal convention is that
    **     the quasi-JD day represents UTC days whether the length is 86399,
    **     86400 or 86401 SI seconds.  In the 1960-1972 era there were
    **     smaller jumps (in either direction) each time the linear UTC(TAI)
    **     expression was changed, and these "mini-leaps" are also included
    **     in the SOFA convention.
    **
    **  5) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale or that are too far in the future
    **     to be trusted.  See iauDat for further details.
    **
    **  6) For calendar conventions and limitations, see iauCal2jd.
    **
    **  Called:
    **     iauJd2cal    JD to Gregorian calendar
    **     iauD2tf      decompose days to hms
    **     iauDat       delta(AT) = TAI-UTC
    **
    **  This revision:  2014 February 15
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    scale = scale.encode()
    iy = c_int()
    im = c_int()
    idd = c_int()
    ihmsf = zeros(shape = (4,1), dtype  = int)
    s = _sofa.iauD2dtf(scale,ndp, d1, d2, byref(iy), byref(im),byref(idd),ihmsf)
    if s != 0:
        warnings.warn(d2dtf_msg[s], UserWarning, 2)
    return iy.value, im.value, idd.value, ihmsf


# iauDat
dat_msg = {0: 'OK', # Unused
                -1:'bad year',
                -2: 'bad month',
                -3: 'bad day(note 3)',
                -4: 'bad fraction(note 4)',
                +1:'dubious year (Note 5)'}
_sofa.iauDat.argtypes = [
    c_int, # iy
    c_int, # im
    c_int, # id
    c_double, # fd
    POINTER(c_double), #delta
]
def dat(iy, im, idd, fd):
    """
    For a given UTC date, calculate Delta(AT) = TAI-UTC.

    Parameters
    ----------
    iy : int
        UTC:  year (Notes 1 and 2).
    im : int
        month (Note 2).
    idd : int
        day (Notes 2 and 3).
    fd : double
        fraction of day (Note 4).

    Returns
    -------
    deltat :  float
        TAI minus UTC, seconds.
    /*
    **  - - - - - - -
    **   i a u D a t
    **  - - - - - -
    **
    **  For a given UTC date, calculate Delta(AT) = TAI-UTC.
    **
    **     :------------------------------------------:
    **     :                                          :
    **     :                 IMPORTANT                :
    **     :                                          :
    **     :  A new version of this function must be  :
    **     :  produced whenever a new leap second is  :
    **     :  announced.  There are four items to     :
    **     :  change on each such occasion:           :
    **     :                                          :
    **     :  1) A new line must be added to the set  :
    **     :     of statements that initialize the    :
    **     :     array "changes".                     :
    **     :                                          :
    **     :  2) The constant IYV must be set to the  :
    **     :     current year.                        :
    **     :                                          :
    **     :  3) The "Latest leap second" comment     :
    **     :     below must be set to the new leap    :
    **     :     second date.                         :
    **     :                                          :
    **     :  4) The "This revision" comment, later,  :
    **     :     must be set to the current date.     :
    **     :                                          :
    **     :  Change (2) must also be carried out     :
    **     :  whenever the function is re-issued,     :
    **     :  even if no leap seconds have been       :
    **     :  added.                                  :
    **     :                                          :
    **     :  Latest leap second:  2016 December 31   :
    **     :                                          :
    **     :__________________________________________:
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  user-replaceable support function.
    **
    **  Given:
    **     iy     int      UTC:  year (Notes 1 and 2)
    **     im     int            month (Note 2)
    **     id     int            day (Notes 2 and 3)
    **     fd     double         fraction of day (Note 4)
    **
    **  Returned:
    **     deltat double   TAI minus UTC, seconds
    **
    **  Returned (function value):
    **            int      status (Note 5):
    **                       1 = dubious year (Note 1)
    **                       0 = OK
    **                      -1 = bad year
    **                      -2 = bad month
    **                      -3 = bad day (Note 3)
    **                      -4 = bad fraction (Note 4)
    **                      -5 = internal error (Note 5)
    **
    **  Notes:
    **
    **  1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
    **     to call the function with an earlier date.  If this is attempted,
    **     zero is returned together with a warning status.
    **
    **     Because leap seconds cannot, in principle, be predicted in
    **     advance, a reliable check for dates beyond the valid range is
    **     impossible.  To guard against gross errors, a year five or more
    **     after the release year of the present function (see the constant
    **     IYV) is considered dubious.  In this case a warning status is
    **     returned but the result is computed in the normal way.
    **
    **     For both too-early and too-late years, the warning status is +1.
    **     This is distinct from the error status -1, which signifies a year
    **     so early that JD could not be computed.
    **
    **  2) If the specified date is for a day which ends with a leap second,
    **     the TAI-UTC value returned is for the period leading up to the
    **     leap second.  If the date is for a day which begins as a leap
    **     second ends, the TAI-UTC returned is for the period following the
    **     leap second.
    **
    **  3) The day number must be in the normal calendar range, for example
    **     1 through 30 for April.  The "almanac" convention of allowing
    **     such dates as January 0 and December 32 is not supported in this
    **     function, in order to avoid confusion near leap seconds.
    **
    **  4) The fraction of day is used only for dates before the
    **     introduction of leap seconds, the first of which occurred at the
    **     end of 1971.  It is tested for validity (0 to 1 is the valid
    **     range) even if not used;  if invalid, zero is used and status -4
    **     is returned.  For many applications, setting fd to zero is
    **     acceptable;  the resulting error is always less than 3 ms (and
    **     occurs only pre-1972).
    **
    **  5) The status value returned in the case where there are multiple
    **     errors refers to the first error detected.  For example, if the
    **     month and day are 13 and 32 respectively, status -2 (bad month)
    **     will be returned.  The "internal error" status refers to a
    **     case that is impossible but causes some compilers to issue a
    **     warning.
    **
    **  6) In cases where a valid result is not available, zero is returned.
    **
    **  References:
    **
    **  1) For dates from 1961 January 1 onwards, the expressions from the
    **     file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
    **
    **  2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
    **     the 1992 Explanatory Supplement.
    **
    **  Called:
    **     iauCal2jd    Gregorian calendar to JD
    **
    **  This revision:  2019 July 5
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    deltat = c_double()
    s = _sofa.iauDat(iy,im,idd,fd,byref(deltat))
    if s != 0:
        warnings.warn(dat_msg[s], UserWarning, 2)
    return deltat.value


# iauDtdb
_sofa.iauDtdb.argtypes = [
    c_double,  # date1
    c_double,  # date2
    c_double,  # ut
    c_double,  # elong
    c_double,  # u
    c_double]  # v
_sofa.iauDtdb.restype =c_double # TDB-TT
def dtdb(date1,date2,ut,elong,u,v):
    """
    An approximation to TDB-TT, the difference between barycentric
**  dynamical time and terrestrial time, for an observer on the Earth.

    Parameters
    ----------
    date1 : float
        date, TDB (Notes 1-3).
    date2 : float
        date, TDB (Notes 1-3).
    ut : float
        universal time (UT1, fraction of one day).
    elong : float
        longitude (east positive, radians).
    u : float
        distance from Earth spin axis (km).
    v : float
        distance north of equatorial plane (km).

    Returns
    -------
    dtdbtt : float
        TDB-TT (seconds).
    /*
    **  - - - - - - - -
    **   i a u D t d b
    **  - - - - - - - -
    **
    **  An approximation to TDB-TT, the difference between barycentric
    **  dynamical time and terrestrial time, for an observer on the Earth.
    **
    **  The different time scales - proper, coordinate and realized - are
    **  related to each other:
    **
    **            TAI             <-  physically realized
    **             :
    **          offset            <-  observed (nominally +32.184s)
    **             :
    **            TT              <-  terrestrial time
    **             :
    **    rate adjustment (L_G)   <-  definition of TT
    **             :
    **            TCG             <-  time scale for GCRS
    **             :
    **      "periodic" terms      <-  iauDtdb  is an implementation
    **             :
    **    rate adjustment (L_C)   <-  function of solar-system ephemeris
    **             :
    **            TCB             <-  time scale for BCRS
    **             :
    **    rate adjustment (-L_B)  <-  definition of TDB
    **             :
    **            TDB             <-  TCB scaled to track TT
    **             :
    **      "periodic" terms      <-  -iauDtdb is an approximation
    **             :
    **            TT              <-  terrestrial time
    **
    **  Adopted values for the various constants can be found in the IERS
    **  Conventions (McCarthy & Petit 2003).
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards Of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     date1,date2   double  date, TDB (Notes 1-3)
    **     ut            double  universal time (UT1, fraction of one day)
    **     elong         double  longitude (east positive, radians)
    **     u             double  distance from Earth spin axis (km)
    **     v             double  distance north of equatorial plane (km)
    **
    **  Returned (function value):
    **                   double  TDB-TT (seconds)
    **
    **  Notes:
    **
    **  1) The date date1+date2 is a Julian Date, apportioned in any
    **     convenient way between the two arguments.  For example,
    **     JD(TT)=2450123.7 could be expressed in any of these ways,
    **     among others:
    **
    **            date1          date2
    **
    **         2450123.7           0.0       (JD method)
    **         2451545.0       -1421.3       (J2000 method)
    **         2400000.5       50123.2       (MJD method)
    **         2450123.5           0.2       (date & time method)
    **
    **     The JD method is the most natural and convenient to use in
    **     cases where the loss of several decimal digits of resolution
    **     is acceptable.  The J2000 method is best matched to the way
    **     the argument is handled internally and will deliver the
    **     optimum resolution.  The MJD method and the date & time methods
    **     are both good compromises between resolution and convenience.
    **
    **     Although the date is, formally, barycentric dynamical time (TDB),
    **     the terrestrial dynamical time (TT) can be used with no practical
    **     effect on the accuracy of the prediction.
    **
    **  2) TT can be regarded as a coordinate time that is realized as an
    **     offset of 32.184s from International Atomic Time, TAI.  TT is a
    **     specific linear transformation of geocentric coordinate time TCG,
    **     which is the time scale for the Geocentric Celestial Reference
    **     System, GCRS.
    **
    **  3) TDB is a coordinate time, and is a specific linear transformation
    **     of barycentric coordinate time TCB, which is the time scale for
    **     the Barycentric Celestial Reference System, BCRS.
    **
    **  4) The difference TCG-TCB depends on the masses and positions of the
    **     bodies of the solar system and the velocity of the Earth.  It is
    **     dominated by a rate difference, the residual being of a periodic
    **     character.  The latter, which is modeled by the present function,
    **     comprises a main (annual) sinusoidal term of amplitude
    **     approximately 0.00166 seconds, plus planetary terms up to about
    **     20 microseconds, and lunar and diurnal terms up to 2 microseconds.
    **     These effects come from the changing transverse Doppler effect
    **     and gravitational red-shift as the observer (on the Earth's
    **     surface) experiences variations in speed (with respect to the
    **     BCRS) and gravitational potential.
    **
    **  5) TDB can be regarded as the same as TCB but with a rate adjustment
    **     to keep it close to TT, which is convenient for many applications.
    **     The history of successive attempts to define TDB is set out in
    **     Resolution 3 adopted by the IAU General Assembly in 2006, which
    **     defines a fixed TDB(TCB) transformation that is consistent with
    **     contemporary solar-system ephemerides.  Future ephemerides will
    **     imply slightly changed transformations between TCG and TCB, which
    **     could introduce a linear drift between TDB and TT;  however, any
    **     such drift is unlikely to exceed 1 nanosecond per century.
    **
    **  6) The geocentric TDB-TT model used in the present function is that of
    **     Fairhead & Bretagnon (1990), in its full form.  It was originally
    **     supplied by Fairhead (private communications with P.T.Wallace,
    **     1990) as a Fortran subroutine.  The present C function contains an
    **     adaptation of the Fairhead code.  The numerical results are
    **     essentially unaffected by the changes, the differences with
    **     respect to the Fairhead & Bretagnon original being at the 1e-20 s
    **     level.
    **
    **     The topocentric part of the model is from Moyer (1981) and
    **     Murray (1983), with fundamental arguments adapted from
    **     Simon et al. 1994.  It is an approximation to the expression
    **     ( v / c ) . ( r / c ), where v is the barycentric velocity of
    **     the Earth, r is the geocentric position of the observer and
    **     c is the speed of light.
    **
    **     By supplying zeroes for u and v, the topocentric part of the
    **     model can be nullified, and the function will return the Fairhead
    **     & Bretagnon result alone.
    **
    **  7) During the interval 1950-2050, the absolute accuracy is better
    **     than +/- 3 nanoseconds relative to time ephemerides obtained by
    **     direct numerical integrations based on the JPL DE405 solar system
    **     ephemeris.
    **
    **  8) It must be stressed that the present function is merely a model,
    **     and that numerical integration of solar-system ephemerides is the
    **     definitive method for predicting the relationship between TCG and
    **     TCB and hence between TT and TDB.
    **
    **  References:
    **
    **     Fairhead, L., & Bretagnon, P., Astron.Astrophys., 229, 240-247
    **     (1990).
    **
    **     IAU 2006 Resolution 3.
    **
    **     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    **     IERS Technical Note No. 32, BKG (2004)
    **
    **     Moyer, T.D., Cel.Mech., 23, 33 (1981).
    **
    **     Murray, C.A., Vectorial Astrometry, Adam Hilger (1983).
    **
    **     Seidelmann, P.K. et al., Explanatory Supplement to the
    **     Astronomical Almanac, Chapter 2, University Science Books (1992).
    **
    **     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    **     Francou, G. & Laskar, J., Astron.Astrophys., 282, 663-683 (1994).
    **
    **  This revision:  2018 January 2
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    dtdbtt  = _sofa.iauDtdb(date1,date2,ut,elong,u,v)
    return dtdbtt


#iauDtf2d
dtf2d_msg = {0: 'OK', # Unused
             3 : 'both of next two' ,
             2 : 'time is after end of day(Note 5)',
             -1:'bad year',
             -2:'bad month',
             -3:'bad day',
             -4:'bad hour',
             -5:'bad minute',
             -6:'bad second(<0)'
             }
try:
    _sofa.iauDtf2d.argtypes = [
        POINTER(c_char), #scale
        c_int, #iy
        c_int, #im
        c_int, # id
        c_int, #ihr
        c_int, #imn
        c_double,#sec
        POINTER(c_double),#d1
        POINTER(c_double),#d2        
        ]
except AttributeError:
    pass
def dtf2d(scale,iy,im,idd,ihr,imn,sec):
    """
Encode date and time fields into 2-part Julian Date (or in the case of UTC a quasi-JD form that includes special provision for leap seconds).
    Parameters
    ----------
    scale : char[]
        time scale ID (Note 1).
    iy : int
        year in Gregorian calendar (Note 2)
    im : int
        month in Gregorian calendar (Note 2)
    idd : int
        day in Gregorian calendar (Note 2).
    ihr : int
        hour.
    imn : int
        minute.
    sec : float
        seconds.

    Returns
    -------
    d1 : float
        double  2-part Julian Date (Notes 3,4).
    d2 : float
        double  2-part Julian Date (Notes 3,4).
    /*
    **  - - - - - - - - -
    **   i a u D t f 2 d
    **  - - - - - - - - -
    **
    **  Encode date and time fields into 2-part Julian Date (or in the case
    **  of UTC a quasi-JD form that includes special provision for leap
    **  seconds).
    **
    **  This function is part of the International Astronomical Union's
    **  SOFA (Standards of Fundamental Astronomy) software collection.
    **
    **  Status:  support function.
    **
    **  Given:
    **     scale     char[]  time scale ID (Note 1)
    **     iy,im,id  int     year, month, day in Gregorian calendar (Note 2)
    **     ihr,imn   int     hour, minute
    **     sec       double  seconds
    **
    **  Returned:
    **     d1,d2     double  2-part Julian Date (Notes 3,4)
    **
    **  Returned (function value):
    **               int     status: +3 = both of next two
    **                               +2 = time is after end of day (Note 5)
    **                               +1 = dubious year (Note 6)
    **                                0 = OK
    **                               -1 = bad year
    **                               -2 = bad month
    **                               -3 = bad day
    **                               -4 = bad hour
    **                               -5 = bad minute
    **                               -6 = bad second (<0)
    **
    **  Notes:
    **
    **  1) scale identifies the time scale.  Only the value "UTC" (in upper
    **     case) is significant, and enables handling of leap seconds (see
    **     Note 4).
    **
    **  2) For calendar conventions and limitations, see iauCal2jd.
    **
    **  3) The sum of the results, d1+d2, is Julian Date, where normally d1
    **     is the Julian Day Number and d2 is the fraction of a day.  In the
    **     case of UTC, where the use of JD is problematical, special
    **     conventions apply:  see the next note.
    **
    **  4) JD cannot unambiguously represent UTC during a leap second unless
    **     special measures are taken.  The SOFA internal convention is that
    **     the quasi-JD day represents UTC days whether the length is 86399,
    **     86400 or 86401 SI seconds.  In the 1960-1972 era there were
    **     smaller jumps (in either direction) each time the linear UTC(TAI)
    **     expression was changed, and these "mini-leaps" are also included
    **     in the SOFA convention.
    **
    **  5) The warning status "time is after end of day" usually means that
    **     the sec argument is greater than 60.0.  However, in a day ending
    **     in a leap second the limit changes to 61.0 (or 59.0 in the case
    **     of a negative leap second).
    **
    **  6) The warning status "dubious year" flags UTCs that predate the
    **     introduction of the time scale or that are too far in the future
    **     to be trusted.  See iauDat for further details.
    **
    **  7) Only in the case of continuous and regular time scales (TAI, TT,
    **     TCG, TCB and TDB) is the result d1+d2 a Julian Date, strictly
    **     speaking.  In the other cases (UT1 and UTC) the result must be
    **     used with circumspection;  in particular the difference between
    **     two such results cannot be interpreted as a precise time
    **     interval.
    **
    **  Called:
    **     iauCal2jd    Gregorian calendar to JD
    **     iauDat       delta(AT) = TAI-UTC
    **     iauJd2cal    JD to Gregorian calendar
    **
    **  This revision:  2013 July 26
    **
    **  SOFA release 2019-07-22
    **
    **  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
    */
    """
    scale = scale.encode()
    d1 = c_double()
    d2 = c_double()
    s = _sofa.iauDtf2d(scale,iy,im,idd,ihr,imn,sec,byref(d1),byref(d2))
    if s != 0:
        warnings.warn(dtf2d_msg[s], UserWarning, 2)
    return d1.value, d2.value


# iauTaitt
try:
    _sofa.iauTaitt.argtypes = [
        c_double,#tai1
        c_double,#tai2
        POINTER(c_double),#tt1
        POINTER(c_double),#tt2       
        ]
except AttributeError:
    pass
def taitt(tai1,tai2):
    """
    Time scale transformation:  International Atomic Time, TAI, to Terrestrial Time, TT.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date.
    tai2 : float
        TAI as a 2-part Julian Date.

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T a i t t
**  - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Terrestrial Time, TT.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tai1+tai2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tai1 is the Julian
**     Day Number and tai2 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tt1 = c_double()
    tt2 = c_double()
    _sofa.iauTaitt(tai1,tai2,byref(tt1),byref(tt2))
    return tt1.value, tt2.value


# iauTaiut1
try:
    _sofa.iauTaiut1.argtypes = [
        c_double,#tai1
        c_double,#tai2
        c_double,#dta
        POINTER(c_double),#ut11
        POINTER(c_double),#ut12    
        ]
except AttributeError:
    pass
def taiut1(tai1,tai2,dta):
    """
    

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date.
    tai2 : float
        TAI as a 2-part Julian Date.
    dta : float
        UT1-TAI in seconds.

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date.
    ut12 : float
        UT1 as a 2-part Julian Date.
/*
**  - - - - - - - - - -
**   i a u T a i u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Universal Time, UT1.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**     dta        double    UT1-TAI in seconds
**
**  Returned:
**     ut11,ut12  double    UT1 as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tai1+tai2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tai1 is the Julian
**     Day Number and tai2 is the fraction of a day.  The returned
**     UT11,UT12 follow suit.
**
**  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
**     available from IERS tabulations.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
**
*/
    """
    ut11 = c_double()
    ut12 = c_double()
    _sofa.iauTaiut1(tai1,tai2,dta,byref(ut11),byref(ut12))
    return ut11.value,ut12.value


# iauTaiutc
Taiutc_msg = {0: 'OK', # Unused
             +1 : 'dubious year (Note 4)' ,
             -1:'unacceptable date',
             }
try:
    _sofa.iauTaiutc.argtypes = [
        c_double,#tai1
        c_double,#tai2
        POINTER(c_double),#utc1
        POINTER(c_double),#utc2    
        ]
except AttributeError:
    pass
def taiutc(tai1,tai2):
    """
    Time scale transformation:  International Atomic Time, TAI, to Coordinated Universal Time, UTC.

    Parameters
    ----------
    tai1 : float
        TAI as a 2-part Julian Date (Note 1)
    tai2 : float
        TAI as a 2-part Julian Date (Note 1)

    Returns
    -------
    utc1 : UTC as a 2-part quasi Julian Date (Notes 1-3)
        DESCRIPTION.
    utc2 : UTC as a 2-part quasi Julian Date (Notes 1-3)
        DESCRIPTION.
/*
**  - - - - - - - - - -
**   i a u T a i u t c
**  - - - - - - - - - -
**
**  Time scale transformation:  International Atomic Time, TAI, to
**  Coordinated Universal Time, UTC.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tai1,tai2  double   TAI as a 2-part Julian Date (Note 1)
**
**  Returned:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-3)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 4)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) tai1+tai2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tai1 is the Julian
**     Day Number and tai2 is the fraction of a day.  The returned utc1
**     and utc2 form an analogous pair, except that a special convention
**     is used, to deal with the problem of leap seconds - see the next
**     note.
**
**  2) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the JD day represents UTC days whether the
**     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
**     there were smaller jumps (in either direction) each time the
**     linear UTC(TAI) expression was changed, and these "mini-leaps"
**     are also included in the SOFA convention.
**
**  3) The function iauD2dtf can be used to transform the UTC quasi-JD
**     into calendar date and clock time, including UTC leap second
**     handling.
**
**  4) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See iauDat for further details.
**
**  Called:
**     iauUtctai    UTC to TAI
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    utc1 = c_double()
    utc2 = c_double()
    s = _sofa.iauTaiutc(tai1,tai2,byref(utc1),byref(utc2))
    if s != 0:
        warnings.warn(Taiutc_msg[s], UserWarning, 2)
    return utc1.value,utc2.value


# iauTcbtdb
try:
    _sofa.iauTcbtdb.argtypes = [
        c_double,#tcb1
        c_double,#tcb2
        POINTER(c_double),#tdb1
        POINTER(c_double),#tdb2    
        ]
except AttributeError:
    pass
def tcbtdb(tcb1,tcb2):
    """
    Time scale transformation:  Barycentric Coordinate Time, TCB, to Barycentric Dynamical Time, TDB.

    Parameters
    ----------
    tcb1 : float
        TCB as a 2-part Julian Date.
    tcb2 : float
        TCB as a 2-part Julian Date.
 Julian Date
**
**  Returned:
**     tdb1,tdb2  double    TDB as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tcb1+tcb2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tcb1 is the Julian
**     Day Number and tcb2 is the fraction of a day.  The returned
**     tdb1,tdb2 follow suit.
**
**  2) The 2006 IAU General Assembly introduced a conventional linear
**     transformation between TDB and TCB.  This transformation
**     compensates for the drift between TCB and terrestrial time TT,
**     and keeps TDB approximately centered on TT.  Because the
**     relationship between TT and TCB depends on the adopted solar
**     system ephemeris, the degree of alignment between TDB and TT over
**     long intervals will vary according to which ephemeris is used.
**     Former definitions of TDB attempted to avoid this problem by
**     stipulating that TDB and TT should differ only by periodic
**     effects.  This is a good description of the nature of the
**     relationship but eluded precise mathematical formulation.  The
**     conventional linear relationship adopted in 2006 sidestepped
**     these difficulties whilst delivering a TDB that in practice was
**     consistent with values before that date.
**
**  3) TDB is essentially the same as Teph, the time argument for the
**     JPL solar system ephemerides.
**
**  Reference:
**
**     IAU 2006 Resolution B3
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    Returns
    -------
    tdb1 : float
        TDB as a 2-part Julian Date.
    tdb2 :float
        TDB as a 2-part Julian Date.

    """
    tdb1 = c_double()
    tdb2 = c_double()
    _sofa.iauTcbtdb(tcb1,tcb2,byref(tdb1),byref(tdb2))
    return tdb1.value,tdb2.value


# iauTcgtt
try:
    _sofa.iauTcgtt.argtypes = [
        c_double,#tcg1
        c_double,#tcg2
        POINTER(c_double),#tt1
        POINTER(c_double),#tt2    
        ]
except AttributeError:
    pass
def tcgtt(tcg1,tcg2):
    """
Time scale transformation:  Geocentric Coordinate Time, TCG, to Terrestrial Time, TT.

    Parameters
    ----------
    tcg1 : float
        TCG as a 2-part Julian Date.
    tcg2 : float
        TCG as a 2-part Julian Date.

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T c g t t
**  - - - - - - - - -
**
**  Time scale transformation:  Geocentric Coordinate Time, TCG, to
**  Terrestrial Time, TT.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tcg1,tcg2  double    TCG as a 2-part Julian Date
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tcg1+tcg2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tcg1 is the Julian
**     Day Number and tcg22 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2000 Resolution B1.9
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tt1 = c_double()
    tt2 = c_double()
    _sofa.iauTcgtt(tcg1,tcg2,byref(tt1),byref(tt2))
    return tt1.value,tt2.value


# iauTdbtcb
try:
    _sofa.iauTdbtcb.argtypes = [
        c_double,#tdb1
        c_double,#tdb2
        POINTER(c_double),#tcb1
        POINTER(c_double),#tcb2 
        ]
except AttributeError:
    pass
def tdbtcb(tdb1,tdb2):
    """
  Time scale transformation:  Barycentric Dynamical Time, TDB, to Barycentric Coordinate Time, TCB.

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date.
    tdb2 : float
        TDB as a 2-part Julian Date.

    Returns
    -------
    tcb1 : float
        TCB as a 2-part Julian Date.
    tcb2 : float
        TCB as a 2-part Julian Date.
/*
**  - - - - - - - - - -
**   i a u T d b t c b
**  - - - - - - - - - -
**
**  Time scale transformation:  Barycentric Dynamical Time, TDB, to
**  Barycentric Coordinate Time, TCB.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tdb1,tdb2  double    TDB as a 2-part Julian Date
**
**  Returned:
**     tcb1,tcb2  double    TCB as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tdb1 is the Julian
**     Day Number and tdb2 is the fraction of a day.  The returned
**     tcb1,tcb2 follow suit.
**
**  2) The 2006 IAU General Assembly introduced a conventional linear
**     transformation between TDB and TCB.  This transformation
**     compensates for the drift between TCB and terrestrial time TT,
**     and keeps TDB approximately centered on TT.  Because the
**     relationship between TT and TCB depends on the adopted solar
**     system ephemeris, the degree of alignment between TDB and TT over
**     long intervals will vary according to which ephemeris is used.
**     Former definitions of TDB attempted to avoid this problem by
**     stipulating that TDB and TT should differ only by periodic
**     effects.  This is a good description of the nature of the
**     relationship but eluded precise mathematical formulation.  The
**     conventional linear relationship adopted in 2006 sidestepped
**     these difficulties whilst delivering a TDB that in practice was
**     consistent with values before that date.
**
**  3) TDB is essentially the same as Teph, the time argument for the
**     JPL solar system ephemerides.
**
**  Reference:
**
**     IAU 2006 Resolution B3
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tcb1 = c_double()
    tcb2 = c_double()
    _sofa.iauTdbtcb(tdb1,tdb2,byref(tcb1),byref(tcb2))
    return tcb1.value,tcb2.value

# iauTdbtt
try:
    _sofa.iauTdbtt.argtypes = [
        c_double,#tdb1
        c_double,#tdb2
        c_double,#dtr
        POINTER(c_double),#tt1
        POINTER(c_double),#tt2 
        ]
except AttributeError:
    pass
def tdbtt(tdb1,tdb2,dtr):
    """
 Time scale transformation:  Barycentric Dynamical Time, TDB, to      Terrestrial Time, TT.   

    Parameters
    ----------
    tdb1 : float
        TDB as a 2-part Julian Date.
    tdb2 : float
        TDB as a 2-part Julian Date.
    dtr : float
        TDB-TT in seconds.

    Returns
    -------
    tt1 : float
       TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T d b t t
**  - - - - - - - - -
**
**  Time scale transformation:  Barycentric Dynamical Time, TDB, to
**  Terrestrial Time, TT.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tdb1,tdb2  double    TDB as a 2-part Julian Date
**     dtr        double    TDB-TT in seconds
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tdb1+tdb2 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where tdb1 is the Julian
**     Day Number and tdb2 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  2) The argument dtr represents the quasi-periodic component of the
**     GR transformation between TT and TCB.  It is dependent upon the
**     adopted solar-system ephemeris, and can be obtained by numerical
**     integration, by interrogating a precomputed time ephemeris or by
**     evaluating a model such as that implemented in the SOFA function
**     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
**     amplitude.
**
**  3) TDB is essentially the same as Teph, the time argument for the
**     JPL solar system ephemerides.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2006 Resolution 3
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
**
*/
    """
    tt1 = c_double()
    tt2 = c_double()
    _sofa.iauTdbtt(tdb1,tdb2,dtr,byref(tt1),byref(tt2))
    return tt1.value,tt2.value


#iauTttai
try:
    _sofa.iauTttai.argtypes = [
        c_double,#tt1
        c_double,#tt2
        POINTER(c_double),#tai1
        POINTER(c_double),#tai2 
        ]
except AttributeError:
    pass
def tttai(tt1,tt2):
    """
Time scale transformation:  Terrestrial Time, TT, to International Atomic Time, TAI.   

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date.
    tai2 : float
        TAI as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T t t a i
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to International
**  Atomic Time, TAI.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned tai1,tai2 follow
**     suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tai1 = c_double()
    tai2 = c_double()
    _sofa.iauTttai(tt1,tt2,byref(tai1),byref(tai2))
    return tai1.value,tai2.value


#iauTttcg
try:
    _sofa.iauTttcg.argtypes = [
        c_double,#tt1
        c_double,#tt2
        POINTER(c_double),#tcg1
        POINTER(c_double),#tcg2 
        ]
except AttributeError:
    pass
def tttcg(tt1,tt2):
    """
   Time scale transformation:  Terrestrial Time, TT, to Geocentric Coordinate Time, TCG. 

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.

    Returns
    -------
    tcg1 : float
        TCG as a 2-part Julian Date.
    tcg2 : float
        TCG as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T t t c g
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Geocentric
**  Coordinate Time, TCG.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned:
**     tcg1,tcg2  double    TCG as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Note:
**
**     tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
**     suit.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2000 Resolution B1.9
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tcg1 = c_double()
    tcg2 = c_double()
    _sofa.iauTttcg(tt1,tt2,byref(tcg1),byref(tcg2))
    return tcg1.value,tcg2.value


# iauTttdb
try:
    _sofa.iauTttdb.argtypes = [
        c_double,#tt1
        c_double,#tt2
        c_double,#dtr
        POINTER(c_double),#tdb1
        POINTER(c_double),#tdb2
        ]
except AttributeError:
    pass
def tttdb(tt1,tt2,dtr):
    """
    Time scale transformation:  Terrestrial Time, TT, to Barycentric Dynamical Time, TDB.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.
    dtr : float
        TDB-TT in seconds.

    Returns
    -------
    tdb1 : float
       TDB as a 2-part Julian Date.
    tdb2 : float
        TDB as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T t t d b
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Barycentric
**  Dynamical Time, TDB.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**     dtr        double    TDB-TT in seconds
**
**  Returned:
**     tdb1,tdb2  double    TDB as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
**     suit.
**
**  2) The argument dtr represents the quasi-periodic component of the
**     GR transformation between TT and TCB.  It is dependent upon the
**     adopted solar-system ephemeris, and can be obtained by numerical
**     integration, by interrogating a precomputed time ephemeris or by
**     evaluating a model such as that implemented in the SOFA function
**     iauDtdb.   The quantity is dominated by an annual term of 1.7 ms
**     amplitude.
**
**  3) TDB is essentially the same as Teph, the time argument for the JPL
**     solar system ephemerides.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     IAU 2006 Resolution 3
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tdb1 = c_double()
    tdb2 = c_double()
    _sofa.iauTttdb(tt1,tt2,dtr,byref(tdb1),byref(tdb2))
    return tdb1.value,tdb2.value


# iauTtut1
try:
    _sofa.iauTtut1.argtypes = [
        c_double,#tt1
        c_double,#tt2
        c_double,#dt
        POINTER(c_double),#ut11
        POINTER(c_double),#ut12
        ]
except AttributeError:
    pass
def ttut1(tt1,tt2,dt):
    """
    Time scale transformation:  Terrestrial Time, TT, to Universal Time, UT1.

    Parameters
    ----------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : float
        TT as a 2-part Julian Date.
    dt : float
        TT-UT1 in seconds.

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date.
    ut12 : float
        UT1 as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u T t u t 1
**  - - - - - - - - -
**
**  Time scale transformation:  Terrestrial Time, TT, to Universal Time,
**  UT1.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     tt1,tt2    double    TT as a 2-part Julian Date
**     dt         double    TT-UT1 in seconds
**
**  Returned:
**     ut11,ut12  double    UT1 as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) tt1+tt2 is Julian Date, apportioned in any convenient way between
**     the two arguments, for example where tt1 is the Julian Day Number
**     and tt2 is the fraction of a day.  The returned ut11,ut12 follow
**     suit.
**
**  2) The argument dt is classical Delta T.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ut11 = c_double()
    ut12 = c_double()
    _sofa.iauTtut1(tt1,tt2,dt,byref(ut11),byref(ut12))
    return ut11.value,ut12.value


#iauUt1tai
try:
    _sofa.iauUt1tai.argtypes = [
        c_double,#ut11
        c_double,#ut12
        c_double,#dta
        POINTER(c_double),#tai1
        POINTER(c_double),#tai2
        ]
except AttributeError:
    pass
def ut1tai(ut11,ut12,dta):
    """
    Time scale transformation:  Universal Time, UT1, to International Atomic Time, TAI.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date.
    ut12 : float
        UT1 as a 2-part Julian Date.
    dta : float
        UT1-TAI in seconds.

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date.
    tai2 : float
        TAI as a 2-part Julian Date.
/*
**  - - - - - - - - - -
**   i a u U t 1 t a i
**  - - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to International
**  Atomic Time, TAI.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     ut11,ut12  double    UT1 as a 2-part Julian Date
**     dta        double    UT1-TAI in seconds
**
**  Returned:
**     tai1,tai2  double    TAI as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) ut11+ut12 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where ut11 is the Julian
**     Day Number and ut12 is the fraction of a day.  The returned
**     tai1,tai2 follow suit.
**
**  2) The argument dta, i.e. UT1-TAI, is an observed quantity, and is
**     available from IERS tabulations.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tai1 = c_double()
    tai2 = c_double()
    _sofa.iauUt1tai(ut11,ut12,dta,byref(tai1),byref(tai2))
    return tai1.value,tai2.value


#iauUt1tt
try:
    _sofa.iauUt1tt.argtypes = [
        c_double,#ut11
        c_double,#ut12
        c_double,#dt
        POINTER(c_double),#tt1
        POINTER(c_double),#tt2
        ]
except AttributeError:
    pass
def ut1tt(ut11,ut12,dt):
    """
    Time scale transformation:  Universal Time, UT1, to Terrestrial Time, TT.

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date.
    ut12 : float
        UT1 as a 2-part Julian Date.
    dt : float
       TT-UT1 in seconds.

    Returns
    -------
    tt1 : float
        TT as a 2-part Julian Date.
    tt2 : foat
       TT as a 2-part Julian Date.
/*
**  - - - - - - - - -
**   i a u U t 1 t t
**  - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to Terrestrial
**  Time, TT.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     ut11,ut12  double    UT1 as a 2-part Julian Date
**     dt         double    TT-UT1 in seconds
**
**  Returned:
**     tt1,tt2    double    TT as a 2-part Julian Date
**
**  Returned (function value):
**                int       status:  0 = OK
**
**  Notes:
**
**  1) ut11+ut12 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where ut11 is the Julian
**     Day Number and ut12 is the fraction of a day.  The returned
**     tt1,tt2 follow suit.
**
**  2) The argument dt is classical Delta T.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    tt1 = c_double()
    tt2 = c_double()
    _sofa.iauUt1tt(ut11,ut12,dt,byref(tt1),byref(tt2))
    return tt1.value,tt2.value


# iauUt1utc
ut1utc_msg = {0: 'OK', # Unused
             +1 : 'dubious year (Note 4)' ,
             -1:'unacceptable date'
             }
try:
    _sofa.iauUt1utc.argtypes = [
        c_double,#ut11
        c_double,#ut12
        c_double,#dut1
        POINTER(c_double),#utc1
        POINTER(c_double),#utc2
        ]
except AttributeError:
    pass
def ut1utc(ut11,ut12,dut1):
    """
    

    Parameters
    ----------
    ut11 : float
        UT1 as a 2-part Julian Date (Note 1).
    ut12 : float
        UT1 as a 2-part Julian Date (Note 1).
    dut1 : float
        Delta UT1: UT1-UTC in seconds (Note 2).

    Returns
    -------
    utc11 : float
        UTC as a 2-part quasi Julian Date (Notes 3,4).
    utc12 : float
        UTC as a 2-part quasi Julian Date (Notes 3,4).
/*
**  - - - - - - - - - -
**   i a u U t 1 u t c
**  - - - - - - - - - -
**
**  Time scale transformation:  Universal Time, UT1, to Coordinated
**  Universal Time, UTC.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 1)
**     dut1       double   Delta UT1: UT1-UTC in seconds (Note 2)
**
**  Returned:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 3,4)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 5)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) ut11+ut12 is Julian Date, apportioned in any convenient way
**     between the two arguments, for example where ut11 is the Julian
**     Day Number and ut12 is the fraction of a day.  The returned utc1
**     and utc2 form an analogous pair, except that a special convention
**     is used, to deal with the problem of leap seconds - see Note 3.
**
**  2) Delta UT1 can be obtained from tabulations provided by the
**     International Earth Rotation and Reference Systems Service.  The
**     value changes abruptly by 1s at a leap second;  however, close to
**     a leap second the algorithm used here is tolerant of the "wrong"
**     choice of value being made.
**
**  3) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the returned quasi JD day UTC1+UTC2 represents
**     UTC days whether the length is 86399, 86400 or 86401 SI seconds.
**
**  4) The function iauD2dtf can be used to transform the UTC quasi-JD
**     into calendar date and clock time, including UTC leap second
**     handling.
**
**  5) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See iauDat for further details.
**
**  Called:
**     iauJd2cal    JD to Gregorian calendar
**     iauDat       delta(AT) = TAI-UTC
**     iauCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    utc1 = c_double()
    utc2 = c_double()
    s = _sofa.iauUt1utc(ut11,ut12,dut1,byref(utc1),byref(utc2))
    if s != 0:
        warnings.warn(ut1utc_msg[s], UserWarning, 2)
    return utc1.value,utc2.value


# iauUtctai
utctai_msg = {0: 'OK', # Unused
             +1 : 'dubious year (Note 4)' ,
             -1:'unacceptable date'
             }
try:
    _sofa.iauUtctai.argtypes = [
        c_double,#utc1
        c_double,#utc2
        POINTER(c_double),#tai1
        POINTER(c_double),#tai2
        ]
except AttributeError:
    pass
def utctai(utc1,utc2):
    """
    Time scale transformation:  Coordinated Universal Time, UTC, to International Atomic Time, TAI

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date (Notes 1-4).
    utc2 : float
        UTC as a 2-part quasi Julian Date (Notes 1-4).

    Returns
    -------
    tai1 : float
        TAI as a 2-part Julian Date (Note 5).
    tai2 : float
        TAI as a 2-part Julian Date (Note 5).
/*
**  - - - - - - - - - -
**   i a u U t c t a i
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  International Atomic Time, TAI.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
**
**  Returned:
**     tai1,tai2  double   TAI as a 2-part Julian Date (Note 5)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 3)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
**     convenient way between the two arguments, for example where utc1
**     is the Julian Day Number and utc2 is the fraction of a day.
**
**  2) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the JD day represents UTC days whether the
**     length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
**     there were smaller jumps (in either direction) each time the
**     linear UTC(TAI) expression was changed, and these "mini-leaps"
**     are also included in the SOFA convention.
**
**  3) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See iauDat for further details.
**
**  4) The function iauDtf2d converts from calendar date and time of day
**     into 2-part Julian Date, and in the case of UTC implements the
**     leap-second-ambiguity convention described above.
**
**  5) The returned TAI1,TAI2 are such that their sum is the TAI Julian
**     Date.
**
**  Called:
**     iauJd2cal    JD to Gregorian calendar
**     iauDat       delta(AT) = TAI-UTC
**     iauCal2jd    Gregorian calendar to JD
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  This revision:  2019 June 20
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
**
*/
    """
    tai1 = c_double()
    tai2 = c_double()
    s = _sofa.iauUtctai(utc1,utc2,byref(tai1),byref(tai2))
    if s != 0:
        warnings.warn(utctai_msg[s], UserWarning, 2)
    return tai1.value,tai2.value


# iauUtcut1
utcut1_msg = {0: 'OK', # Unused
             +1 : 'dubious year (Note 4)' ,
             -1:'unacceptable date'
             }
try:
    _sofa.iauUtcut1.argtypes = [
        c_double,#utc1
        c_double,#utc2
        c_double,#dut1
        POINTER(c_double),#ut11
        POINTER(c_double),#ut12
        ]
except AttributeError:
    pass
def utcut1(utc1,utc2,dut1):
    """
    Time scale transformation:  Coordinated Universal Time, UTC, to Universal Time, UT1.

    Parameters
    ----------
    utc1 : float
        UTC as a 2-part quasi Julian Date (Notes 1-4).
    utc2 : float
        UTC as a 2-part quasi Julian Date (Notes 1-4).
    dut1 : float
        Delta UT1 = UT1-UTC in seconds (Note 5).

    Returns
    -------
    ut11 : float
        UT1 as a 2-part Julian Date (Note 6).
    ut12 : float
        UT1 as a 2-part Julian Date (Note 6).
/*
**  - - - - - - - - - -
**   i a u U t c u t 1
**  - - - - - - - - - -
**
**  Time scale transformation:  Coordinated Universal Time, UTC, to
**  Universal Time, UT1.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     utc1,utc2  double   UTC as a 2-part quasi Julian Date (Notes 1-4)
**     dut1       double   Delta UT1 = UT1-UTC in seconds (Note 5)
**
**  Returned:
**     ut11,ut12  double   UT1 as a 2-part Julian Date (Note 6)
**
**  Returned (function value):
**                int      status: +1 = dubious year (Note 3)
**                                  0 = OK
**                                 -1 = unacceptable date
**
**  Notes:
**
**  1) utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
**     convenient way between the two arguments, for example where utc1
**     is the Julian Day Number and utc2 is the fraction of a day.
**
**  2) JD cannot unambiguously represent UTC during a leap second unless
**     special measures are taken.  The convention in the present
**     function is that the JD day represents UTC days whether the
**     length is 86399, 86400 or 86401 SI seconds.
**
**  3) The warning status "dubious year" flags UTCs that predate the
**     introduction of the time scale or that are too far in the future
**     to be trusted.  See iauDat for further details.
**
**  4) The function iauDtf2d converts from calendar date and time of
**     day into 2-part Julian Date, and in the case of UTC implements
**     the leap-second-ambiguity convention described above.
**
**  5) Delta UT1 can be obtained from tabulations provided by the
**     International Earth Rotation and Reference Systems Service.
**     It is the caller's responsibility to supply a dut1 argument
**     containing the UT1-UTC value that matches the given UTC.
**
**  6) The returned ut11,ut12 are such that their sum is the UT1 Julian
**     Date.
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**  Called:
**     iauJd2cal    JD to Gregorian calendar
**     iauDat       delta(AT) = TAI-UTC
**     iauUtctai    UTC to TAI
**     iauTaiut1    TAI to UT1
**
**  This revision:  2013 August 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ut11 = c_double()
    ut12 = c_double()
    s = _sofa.iauUtcut1(utc1,utc2,dut1,byref(ut11),byref(ut12))
    if s != 0:
        warnings.warn(utcut1_msg[s], UserWarning, 2)
    return ut11.value,ut12.value

###############################################################################################################################
#####################################################################        Earth rotation angle and sidereal time
################################################################

# iauEe00
try:
    _sofa.iauEe00.argtypes = [
        c_double,#date1
        c_double,#date2
        c_double,#epsa
        c_double,#dpsi
        ]
    _sofa.iauEe00.restype = c_double
except AttributeError:
    pass
def ee00(date1,date2,epsa,dpsi):
    """
    The equation of the equinoxes, compatible with IAU 2000 resolutions,given the nutation in longitude and the mean obliquity.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    epsa : float
        mean obliquity (Note 2).
    dpsi : float
        nutation in longitude (Note 3).

    Returns
    -------
    ee : float
        equation of the equinoxes (Note 4).
/*
**  - - - - - - - -
**   i a u E e 0 0
**  - - - - - - - -
**
**  The equation of the equinoxes, compatible with IAU 2000 resolutions,
**  given the nutation in longitude and the mean obliquity.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**     epsa         double    mean obliquity (Note 2)
**     dpsi         double    nutation in longitude (Note 3)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The obliquity, in radians, is mean of date.
**
**  3) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  4) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     iauEect00    equation of the equinoxes complementary terms
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2008 May 16
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ee = _sofa.iauEe00(date1,date2,epsa,dpsi)
    return ee


# iauEe00a
try:
    _sofa.iauEe00a.argtypes = [
        c_double,#date1
        c_double,#date2
        ]
    _sofa.iauEe00a.restype = c_double
except AttributeError:
    pass
def ee00a(date1,date2):
    """
    Equation of the equinoxes, compatible with IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    ee : float
        equation of the equinoxes (Note 2).
/*
**  - - - - - - - - -
**   i a u E e 0 0 a
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  3) The result is compatible with the IAU 2000 resolutions.  For
**     further details, see IERS Conventions 2003 and Capitaine et al.
**     (2002).
**
**  Called:
**     iauPr00      IAU 2000 precession adjustments
**     iauObl80     mean obliquity, IAU 1980
**     iauNut00a    nutation, IAU 2000A
**     iauEe00      equation of the equinoxes, IAU 2000
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003).
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004).
**
**  This revision:  2008 May 16
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ee = _sofa.iauEe00a(date1,date2)
    return ee


# iauEe00b
try:
    _sofa.iauEe00b.argtypes = [
        c_double,#date1
        c_double,#date2
        ]
    _sofa.iauEe00b.restype = c_double
except AttributeError:
    pass
def ee00b(date1,date2):
    """
    Equation of the equinoxes, compatible with IAU 2000 resolutions but using the truncated nutation model IAU 2000B.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    ee : float
        equation of the equinoxes (Note 2).
/*
**  - - - - - - - - -
**   i a u E e 0 0 b
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions but
**  using the truncated nutation model IAU 2000B.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  3) The result is compatible with the IAU 2000 resolutions except
**     that accuracy has been compromised for the sake of speed.  For
**     further details, see McCarthy & Luzum (2001), IERS Conventions
**     2003 and Capitaine et al. (2003).
**
**  Called:
**     iauPr00      IAU 2000 precession adjustments
**     iauObl80     mean obliquity, IAU 1980
**     iauNut00b    nutation, IAU 2000B
**     iauEe00      equation of the equinoxes, IAU 2000
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
**     precession-nutation of the celestial pole", Celestial Mechanics &
**     Dynamical Astronomy, 85, 37-49 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2008 May 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ee = _sofa.iauEe00b(date1,date2)
    return ee


# iauEe06a
try:
    _sofa.iauEe06a.argtypes = [
        c_double,#date1
        c_double,#date2
        ]
    _sofa.iauEe06a.restype = c_double
except AttributeError:
    pass
def ee06a(date1,date2):
    """
    Equation of the equinoxes, compatible with IAU 2000 resolutions and IAU 2006/2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    ee : float
        equation of the equinoxes (Note 2).
/*
**  - - - - - - - - -
**   i a u E e 0 6 a
**  - - - - - - - - -
**
**  Equation of the equinoxes, compatible with IAU 2000 resolutions and
**  IAU 2006/2000A precession-nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  Called:
**     iauAnpm      normalize angle into range +/- pi
**     iauGst06a    Greenwich apparent sidereal time, IAU 2006/2000A
**     iauGmst06    Greenwich mean sidereal time, IAU 2006
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2008 May 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ee = _sofa.iauEe06a(date1,date2)
    return ee


# iauEect00
try:
    _sofa.iauEect00.argtypes = [
        c_double,#date1
        c_double,#date2
        ]
    _sofa.iauEect00.restype = c_double
except AttributeError:
    pass
def eect00(date1,date2):
    """
    Equation of the equinoxes complementary terms, consistent with IAU 2000 resolutions.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    eect : float
        complementary terms (Note 2).
/*
**  - - - - - - - - - -
**   i a u E e c t 0 0
**  - - - - - - - - - -
**
**  Equation of the equinoxes complementary terms, consistent with
**  IAU 2000 resolutions.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   complementary terms (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The "complementary terms" are part of the equation of the
**     equinoxes (EE), classically the difference between apparent and
**     mean Sidereal Time:
**
**        GAST = GMST + EE
**
**     with:
**
**        EE = dpsi * cos(eps)
**
**     where dpsi is the nutation in longitude and eps is the obliquity
**     of date.  However, if the rotation of the Earth were constant in
**     an inertial frame the classical formulation would lead to
**     apparent irregularities in the UT1 timescale traceable to side-
**     effects of precession-nutation.  In order to eliminate these
**     effects from UT1, "complementary terms" were introduced in 1994
**     (IAU, 1994) and took effect from 1997 (Capitaine and Gontier,
**     1993):
**
**        GAST = GMST + CT + EE
**
**     By convention, the complementary terms are included as part of
**     the equation of the equinoxes rather than as part of the mean
**     Sidereal Time.  This slightly compromises the "geometrical"
**     interpretation of mean sidereal time but is otherwise
**     inconsequential.
**
**     The present function computes CT in the above expression,
**     compatible with IAU 2000 resolutions (Capitaine et al., 2002, and
**     IERS Conventions 2003).
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFalp03    mean anomaly of the Sun
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFad03     mean elongation of the Moon from the Sun
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N. & Gontier, A.-M., Astron.Astrophys., 275,
**     645-650 (1993)
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astron.Astrophys., 406,
**     1135-1149 (2003)
**
**     IAU Resolution C7, Recommendation 3 (1994)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eect = _sofa.iauEect00(date1,date2)
    return eect


# iauEqeq94
try:
    _sofa.iauEqeq94.argtypes = [
        c_double,#date1
        c_double,#date2
        ]
    _sofa.iauEqeq94.restype = c_double
except AttributeError:
    pass
def eqeq94(date1,date2):
    """
    Equation of the equinoxes, IAU 1994 model.

    Parameters
    ----------
    date1 : float
        TDB date (Note 1).
    date2 : float
        TDB date (Note 1).

    Returns
    -------
    ee : float
        equation of the equinoxes (Note 2).
/*
**  - - - - - - - - - -
**   i a u E q e q 9 4
**  - - - - - - - - - -
**
**  Equation of the equinoxes, IAU 1994 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double     TDB date (Note 1)
**
**  Returned (function value):
**                   double     equation of the equinoxes (Note 2)
**
**  Notes:
**
**  1) The date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result, which is in radians, operates in the following sense:
**
**        Greenwich apparent ST = GMST + equation of the equinoxes
**
**  Called:
**     iauAnpm      normalize angle into range +/- pi
**     iauNut80     nutation, IAU 1980
**     iauObl80     mean obliquity, IAU 1980
**
**  References:
**
**     IAU Resolution C7, Recommendation 3 (1994).
**
**     Capitaine, N. & Gontier, A.-M., 1993, Astron.Astrophys., 275,
**     645-650.
**
**  This revision:  2017 October 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ee = _sofa.iauEqeq94(date1,date2)
    return ee


# iauEra00
try:
    _sofa.iauEra00.argtypes = [
        c_double,#dj1
        c_double,#dj2
        ]
    _sofa.iauEra00.restype = c_double
except AttributeError:
    pass
def era00(jd1,jd2):
    """
    Earth rotation angle (IAU 2000 model).

    Parameters
    ----------
    jd1 : float
        UT1 as a 2-part Julian Date (see note)
    jd2 : float
        UT1 as a 2-part Julian Date (see note)

    Returns
    -------
    theta : float
        Earth rotation angle (radians), range 0-2pi.
/*
**  - - - - - - - - -
**   i a u E r a 0 0
**  - - - - - - - - -
**
**  Earth rotation angle (IAU 2000 model).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     dj1,dj2   double    UT1 as a 2-part Julian Date (see note)
**
**  Returned (function value):
**               double    Earth rotation angle (radians), range 0-2pi
**
**  Notes:
**
**  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
**     convenient way between the arguments dj1 and dj2.  For example,
**     JD(UT1)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             dj1            dj2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  The date & time method is
**     best matched to the algorithm used:  maximum precision is
**     delivered when the dj1 argument is for 0hrs UT1 on the day in
**     question and the dj2 argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) The algorithm is adapted from Expression 22 of Capitaine et al.
**     2000.  The time argument has been expressed in days directly,
**     and, to retain precision, integer contributions have been
**     eliminated.  The same formulation is given in IERS Conventions
**     (2003), Chap. 5, Eq. 14.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine N., Guinot B. and McCarthy D.D, 2000, Astron.
**     Astrophys., 355, 398-405.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    theta = _sofa.iauEra00(jd1,jd2)
    return theta


# iauGmst00
try:
    _sofa.iauGmst00.argtypes = [
        c_double,#uta
        c_double,#utb
        c_double,#tta
        c_double,#ttb
        ]
    _sofa.iauGmst00.restype = c_double
except AttributeError:
    pass
def gmst00(uta,utb,tta,ttb):
    """
    Greenwich mean sidereal time (model consistent with IAU 2000 resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    tta : float
        TT as a 2-part Julian Date (Notes 1,2).
    ttb : float
        TT as a 2-part Julian Date (Notes 1,2).

    Returns
    -------
    gmst : float
        Greenwich mean sidereal time (radians).
/*
**  - - - - - - - - - -
**   i a u G m s t 0 0
**  - - - - - - - - - -
**
**  Greenwich mean sidereal time (model consistent with IAU 2000
**  resolutions).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich mean sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A         Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession.  If UT1 is used for
**     both purposes, errors of order 100 microarcseconds result.
**
**  3) This GMST is compatible with the IAU 2000 resolutions and must be
**     used only in conjunction with other IAU 2000 compatible
**     components such as precession-nutation and equation of the
**     equinoxes.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  5) The algorithm is from Capitaine et al. (2003) and IERS
**     Conventions 2003.
**
**  Called:
**     iauEra00     Earth rotation angle, IAU 2000
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gmst = _sofa.iauGmst00(uta,utb,tta,ttb)
    return gmst


# iauGmst06
try:
    _sofa.iauGmst06.argtypes = [
        c_double,#uta
        c_double,#utb
        c_double,#tta
        c_double,#ttb
        ]
    _sofa.iauGmst06.restype = c_double
except AttributeError:
    pass
def gmst06(uta,utb,tta,ttb):
    """
    Greenwich mean sidereal time (consistent with IAU 2006 precession).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    tta : float
        TT as a 2-part Julian Date (Notes 1,2).
    ttb : float
        TT as a 2-part Julian Date (Notes 1,2).

    Returns
    -------
    gmst : float
        Greenwich mean sidereal time (radians).
    /*
**  - - - - - - - - - -
**   i a u G m s t 0 6
**  - - - - - - - - - -
**
**  Greenwich mean sidereal time (consistent with IAU 2006 precession).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich mean sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A        Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     rotation angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession.  If UT1 is used for
**     both purposes, errors of order 100 microarcseconds result.
**
**  3) This GMST is compatible with the IAU 2006 precession and must not
**     be used with other precession models.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  Called:
**     iauEra00     Earth rotation angle, IAU 2000
**     iauAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2005,
**     Astron.Astrophys. 432, 355
**
**  This revision:  2013 June 18
** 
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/

    """
    gmst = _sofa.iauGmst06(uta,utb,tta,ttb)
    return gmst


# iauGmst82
try:
    _sofa.iauGmst82.argtypes = [
        c_double,#dj1
        c_double,#dj2
        ]
    _sofa.iauGmst82.restype = c_double
except AttributeError:
    pass
def gmst82(dj1,dj2):
    """
    Universal Time to Greenwich mean sidereal time (IAU 1982 model).

    Parameters
    ----------
    dj1 : float
        UT1 Julian Date (see note).
    dj2 : float
        UT1 Julian Date (see note).

    Returns
    -------
    gmst : float
        Greenwich mean sidereal time (radians).
/*
**  - - - - - - - - - -
**   i a u G m s t 8 2
**  - - - - - - - - - -
**
**  Universal Time to Greenwich mean sidereal time (IAU 1982 model).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     dj1,dj2    double    UT1 Julian Date (see note)
**
**  Returned (function value):
**                double    Greenwich mean sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 date dj1+dj2 is a Julian Date, apportioned in any
**     convenient way between the arguments dj1 and dj2.  For example,
**     JD(UT1)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             dj1            dj2
**
**         2450123.7          0          (JD method)
**          2451545        -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5         0.2         (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  The date & time method is
**     best matched to the algorithm used:  maximum accuracy (or, at
**     least, minimum noise) is delivered when the dj1 argument is for
**     0hrs UT1 on the day in question and the dj2 argument lies in the
**     range 0 to 1, or vice versa.
**
**  2) The algorithm is based on the IAU 1982 expression.  This is
**     always described as giving the GMST at 0 hours UT1.  In fact, it
**     gives the difference between the GMST and the UT, the steady
**     4-minutes-per-day drawing-ahead of ST with respect to UT.  When
**     whole days are ignored, the expression happens to equal the GMST
**     at 0 hours UT1 each day.
**
**  3) In this function, the entire UT1 (the sum of the two arguments
**     dj1 and dj2) is used directly as the argument for the standard
**     formula, the constant term of which is adjusted by 12 hours to
**     take account of the noon phasing of Julian Date.  The UT1 is then
**     added, but omitting whole days to conserve accuracy.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Transactions of the International Astronomical Union,
**     XVIII B, 67 (1983).
**
**     Aoki et al., Astron.Astrophys., 105, 359-361 (1982).
**
**  This revision:  2017 October 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gmst = _sofa.iauGmst82(dj1,dj2)
    return gmst


# iauGst00a
try:
    _sofa.iauGst00a.argtypes = [
        c_double,#uta
        c_double,#utb
        c_double,#tta
        c_double,#ttb
        ]
    _sofa.iauGst00a.restype = c_double
except AttributeError:
    pass
def gst00a(uta,utb,tta,ttb):
    """
    Greenwich apparent sidereal time (consistent with IAU 2000 resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    tta : float
        TT as a 2-part Julian Date (Notes 1,2).
    ttb : float
        TT as a 2-part Julian Date (Notes 1,2).

    Returns
    -------
    gst : float
        Greenwich apparent sidereal time (radians).
/*
**  - - - - - - - - - -
**   i a u G s t 0 0 a
**  - - - - - - - - - -
**
**  Greenwich apparent sidereal time (consistent with IAU 2000
**  resolutions).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A        Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession-nutation.  If UT1 is
**     used for both purposes, errors of order 100 microarcseconds
**     result.
**
**  3) This GAST is compatible with the IAU 2000 resolutions and must be
**     used only in conjunction with other IAU 2000 compatible
**     components such as precession-nutation.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  5) The algorithm is from Capitaine et al. (2003) and IERS
**     Conventions 2003.
**
**  Called:
**     iauGmst00    Greenwich mean sidereal time, IAU 2000
**     iauEe00a     equation of the equinoxes, IAU 2000A
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gst = _sofa.iauGst00a(uta,utb,tta,ttb)
    return gst


# iauGst00b
try:
    _sofa.iauGst00b.argtypes = [
        c_double,#uta
        c_double,#utb
        ]
    _sofa.iauGst00b.restype = c_double
except AttributeError:
    pass
def gst00b(uta,utb):
    """
    

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).

    Returns
    -------
    gst : float
        Greenwich apparent sidereal time (radians).
/*
**  - - - - - - - - - -
**   i a u G s t 0 0 b
**  - - - - - - - - - -
**
**  Greenwich apparent sidereal time (consistent with IAU 2000
**  resolutions but using the truncated nutation model IAU 2000B).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 date uta+utb is a Julian Date, apportioned in any
**     convenient way between the argument pair.  For example,
**     JD=2450123.7 could be expressed in any of these ways, among
**     others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) The result is compatible with the IAU 2000 resolutions, except
**     that accuracy has been compromised for the sake of speed and
**     convenience in two respects:
**
**     . UT is used instead of TDB (or TT) to compute the precession
**       component of GMST and the equation of the equinoxes.  This
**       results in errors of order 0.1 mas at present.
**
**     . The IAU 2000B abridged nutation model (McCarthy & Luzum, 2001)
**       is used, introducing errors of up to 1 mas.
**
**  3) This GAST is compatible with the IAU 2000 resolutions and must be
**     used only in conjunction with other IAU 2000 compatible
**     components such as precession-nutation.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  5) The algorithm is from Capitaine et al. (2003) and IERS
**     Conventions 2003.
**
**  Called:
**     iauGmst00    Greenwich mean sidereal time, IAU 2000
**     iauEe00b     equation of the equinoxes, IAU 2000B
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Capitaine, N., Wallace, P.T. and McCarthy, D.D., "Expressions to
**     implement the IAU 2000 definition of UT1", Astronomy &
**     Astrophysics, 406, 1135-1149 (2003)
**
**     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
**     precession-nutation of the celestial pole", Celestial Mechanics &
**     Dynamical Astronomy, 85, 37-49 (2003)
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gst = _sofa.iauGst00b(uta,utb)
    return gst


# iauGst06
try:
    _sofa.iauGst06.argtypes = [
        c_double,#uta
        c_double,#utb
        c_double,#tta
        c_double,#ttb
        ndpointer(shape=(3,3), dtype = c_double)
        ]
    _sofa.iauGst06.restype = c_double
except AttributeError:
    pass
def gst06(uta,utb,tta,ttb,rnpb):
    """
    
Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    tta : float
        TT as a 2-part Julian Date (Notes 1,2).
    ttb : float
        TT as a 2-part Julian Date (Notes 1,2).
    rnpb : float
        nutation x precession x bias matrix.

    Returns
    -------
    gst : float
        Greenwich apparent sidereal time (radians).
/*
**  - - - - - - - - -
**   i a u G s t 0 6
**  - - - - - - - - -
**
**  Greenwich apparent sidereal time, IAU 2006, given the NPB matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     uta,utb  double        UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb  double        TT as a 2-part Julian Date (Notes 1,2)
**     rnpb     double[3][3]  nutation x precession x bias matrix
**
**  Returned (function value):
**              double        Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A        Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     rotation angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession-nutation.  If UT1 is
**     used for both purposes, errors of order 100 microarcseconds
**     result.
**
**  3) Although the function uses the IAU 2006 series for s+XY/2, it is
**     otherwise independent of the precession-nutation model and can in
**     practice be used with any equinox-based NPB matrix.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  Called:
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS06       the CIO locator s, given X,Y, IAU 2006
**     iauAnp       normalize angle into range 0 to 2pi
**     iauEra00     Earth rotation angle, IAU 2000
**     iauEors      equation of the origins, given NPB matrix and s
**
**  Reference:
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gst = _sofa.iauGst06(uta,utb,tta,ttb,rnpb)
    return gst


# iauGst06a
try:
    _sofa.iauGst06a.argtypes = [
        c_double,#uta
        c_double,#utb
        c_double,#tta
        c_double,#ttb
        ]
    _sofa.iauGst06a.restype = c_double
except AttributeError:
    pass
def gst06a(uta,utb,tta,ttb):
    """
    Greenwich apparent sidereal time (consistent with IAU 2000 and 2006 resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    tta : float
       TT as a 2-part Julian Date (Notes 1,2).
    ttb : float
        TT as a 2-part Julian Date (Notes 1,2)

    Returns
    -------
    gst : float
        Greenwich apparent sidereal time (radians).
/*
**  - - - - - - - - - -
**   i a u G s t 0 6 a
**  - - - - - - - - - -
**
**  Greenwich apparent sidereal time (consistent with IAU 2000 and 2006
**  resolutions).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**     tta,ttb    double    TT as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 and TT dates uta+utb and tta+ttb respectively, are both
**     Julian Dates, apportioned in any convenient way between the
**     argument pairs.  For example, JD=2450123.7 could be expressed in
**     any of these ways, among others:
**
**            Part A        Part B
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable (in the case of UT;  the TT is not at all critical
**     in this respect).  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     rotation angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) Both UT1 and TT are required, UT1 to predict the Earth rotation
**     and TT to predict the effects of precession-nutation.  If UT1 is
**     used for both purposes, errors of order 100 microarcseconds
**     result.
**
**  3) This GAST is compatible with the IAU 2000/2006 resolutions and
**     must be used only in conjunction with IAU 2006 precession and
**     IAU 2000A nutation.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  Called:
**     iauPnm06a    classical NPB matrix, IAU 2006/2000A
**     iauGst06     Greenwich apparent ST, IAU 2006, given NPB matrix
**
**  Reference:
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gst = _sofa.iauGst06a(uta,utb,tta,ttb)
    return gst


# iauGst94
try:
    _sofa.iauGst94.argtypes = [
        c_double,#uta
        c_double,#utb
        ]
    _sofa.iauGst94.restype = c_double
except AttributeError:
    pass
def gst94(uta,utb):
    """
    Greenwich apparent sidereal time (consistent with IAU 1982/94 resolutions).

    Parameters
    ----------
    uta : float
        UT1 as a 2-part Julian Date (Notes 1,2).
    utb : float
        UT1 as a 2-part Julian Date (Notes 1,2).

    Returns
    -------
    gst : float
        Greenwich apparent sidereal time (radians).
/*
**  - - - - - - - - -
**   i a u G s t 9 4
**  - - - - - - - - -
**
**  Greenwich apparent sidereal time (consistent with IAU 1982/94
**  resolutions).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     uta,utb    double    UT1 as a 2-part Julian Date (Notes 1,2)
**
**  Returned (function value):
**                double    Greenwich apparent sidereal time (radians)
**
**  Notes:
**
**  1) The UT1 date uta+utb is a Julian Date, apportioned in any
**     convenient way between the argument pair.  For example,
**     JD=2450123.7 could be expressed in any of these ways, among
**     others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  For UT, the date & time
**     method is best matched to the algorithm that is used by the Earth
**     Rotation Angle function, called internally:  maximum precision is
**     delivered when the uta argument is for 0hrs UT1 on the day in
**     question and the utb argument lies in the range 0 to 1, or vice
**     versa.
**
**  2) The result is compatible with the IAU 1982 and 1994 resolutions,
**     except that accuracy has been compromised for the sake of
**     convenience in that UT is used instead of TDB (or TT) to compute
**     the equation of the equinoxes.
**
**  3) This GAST must be used only in conjunction with contemporaneous
**     IAU standards such as 1976 precession, 1980 obliquity and 1982
**     nutation.  It is not compatible with the IAU 2000 resolutions.
**
**  4) The result is returned in the range 0 to 2pi.
**
**  Called:
**     iauGmst82    Greenwich mean sidereal time, IAU 1982
**     iauEqeq94    equation of the equinoxes, IAU 1994
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992)
**
**     IAU Resolution C7, Recommendation 3 (1994)
**
**  This revision:  2008 May 16
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gst = _sofa.iauGst94(uta,utb)
    return gst


############################################################################  Ephemerides
################################################################

#iauEpv00
epv00_msg = {0: 'OK', # Unused
             +1 : 'warning: date outside the range 1900-2100 AD'}
try:
    _sofa.iauEpv00.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(2,3),dtype=float),#pvh
        ndpointer(shape=(2,3),dtype=float)]#pvb
except AttributeError:
    pass
def epv00(date1,date2):
    """
    Earth position and velocity, heliocentric and barycentric, with  respect to the Barycentric Celestial Reference System.

    Parameters
    ----------
    date1 : float
        TDB date (Note 1).
    date2 : float
        TDB date (Note 1).

    Returns
    -------
    pvh : double[2][3]
        heliocentric Earth position/velocity.
    pvb : double[2][3]
        barycentric Earth position/velocity.
/*
**  - - - - - - - - -
**   i a u E p v 0 0
**  - - - - - - - - -
**
**  Earth position and velocity, heliocentric and barycentric, with
**  respect to the Barycentric Celestial Reference System.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double        TDB date (Note 1)
**
**  Returned:
**     pvh          double[2][3]  heliocentric Earth position/velocity
**     pvb          double[2][3]  barycentric Earth position/velocity
**
**  Returned (function value):
**                  int           status: 0 = OK
**                                       +1 = warning: date outside
**                                            the range 1900-2100 AD
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways, among
**     others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  However,
**     the accuracy of the result is more likely to be limited by the
**     algorithm itself than the way the date has been expressed.
**
**     n.b. TT can be used instead of TDB in most applications.
**
**  2) On return, the arrays pvh and pvb contain the following:
**
**        pvh[0][0]  x       }
**        pvh[0][1]  y       } heliocentric position, au
**        pvh[0][2]  z       }
**
**        pvh[1][0]  xdot    }
**        pvh[1][1]  ydot    } heliocentric velocity, au/d
**        pvh[1][2]  zdot    }
**
**        pvb[0][0]  x       }
**        pvb[0][1]  y       } barycentric position, au
**        pvb[0][2]  z       }
**
**        pvb[1][0]  xdot    }
**        pvb[1][1]  ydot    } barycentric velocity, au/d
**        pvb[1][2]  zdot    }
**
**     The vectors are with respect to the Barycentric Celestial
**     Reference System.  The time unit is one day in TDB.
**
**  3) The function is a SIMPLIFIED SOLUTION from the planetary theory
**     VSOP2000 (X. Moisson, P. Bretagnon, 2001, Celes. Mechanics &
**     Dyn. Astron., 80, 3/4, 205-213) and is an adaptation of original
**     Fortran code supplied by P. Bretagnon (private comm., 2000).
**
**  4) Comparisons over the time span 1900-2100 with this simplified
**     solution and the JPL DE405 ephemeris give the following results:
**
**                                RMS    max
**           Heliocentric:
**              position error    3.7   11.2   km
**              velocity error    1.4    5.0   mm/s
**
**           Barycentric:
**              position error    4.6   13.4   km
**              velocity error    1.4    4.9   mm/s
**
**     Comparisons with the JPL DE406 ephemeris show that by 1800 and
**     2200 the position errors are approximately double their 1900-2100
**     size.  By 1500 and 2500 the deterioration is a factor of 10 and
**     by 1000 and 3000 a factor of 60.  The velocity accuracy falls off
**     at about half that rate.
**
**  5) It is permissible to use the same array for pvh and pvb, which
**     will receive the barycentric values.
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    pvh = zeros(shape = (2,3), dtype  = float)
    pvb = zeros(shape = (2,3), dtype  = float)
    s = _sofa.iauEpv00(date1,date2,pvh,pvb)
    if s != 0:
        warnings.warn(epv00_msg[s], UserWarning, 2)
    return pvh, pvb


#iauPlan94
plan94_msg = {0: 'OK', # Unused
             -1 : 'illegal NP (outside 1-8)',
             +1 : 'warning: year outside 1000-3000',
             +2 : 'warning: failed to converge'}
try:
    _sofa.iauPlan94.argtypes = [
        c_double,#date1
        c_double,#date2
        c_int,#np
        ndpointer(shape=(2,3),dtype=float)]#pv
except AttributeError:
    pass
def plan94(date1,date2,np):
    """
  This function is part of the International Astronomical Union's SOFA (Standards Of Fundamental Astronomy) software collection.  

    Parameters
    ----------
    date1 : float
        TDB date part A (Note 1).
    date2 : float
        TDB date part B (Note 1).
    np : int
        planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,5=Jupiter, 6=Saturn, 7=Uranus,8=Neptune).

    Returns
    -------
    pv : float
        double[2][3] planet p,v (heliocentric, J2000.0, au,au/d).
/*
**  - - - - - - - - - -
**   i a u P l a n 9 4
**  - - - - - - - - - -
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Approximate heliocentric position and velocity of a nominated major
**  planet:  Mercury, Venus, EMB, Mars, Jupiter, Saturn, Uranus or
**  Neptune (but not the Earth itself).
**
**  Given:
**     date1  double       TDB date part A (Note 1)
**     date2  double       TDB date part B (Note 1)
**     np     int          planet (1=Mercury, 2=Venus, 3=EMB, 4=Mars,
**                             5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune)
**
**  Returned (argument):
**     pv     double[2][3] planet p,v (heliocentric, J2000.0, au,au/d)
**
**  Returned (function value):
**            int          status: -1 = illegal NP (outside 1-8)
**                                  0 = OK
**                                 +1 = warning: year outside 1000-3000
**                                 +2 = warning: failed to converge
**
**  Notes:
**
**  1) The date date1+date2 is in the TDB time scale (in practice TT can
**     be used) and is a Julian Date, apportioned in any convenient way
**     between the two arguments.  For example, JD(TDB)=2450123.7 could
**     be expressed in any of these ways, among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     resolution.  The MJD method and the date & time methods are both
**     good compromises between resolution and convenience.  The limited
**     accuracy of the present algorithm is such that any of the methods
**     is satisfactory.
**
**  2) If an np value outside the range 1-8 is supplied, an error status
**     (function value -1) is returned and the pv vector set to zeroes.
**
**  3) For np=3 the result is for the Earth-Moon Barycenter.  To obtain
**     the heliocentric position and velocity of the Earth, use instead
**     the SOFA function iauEpv00.
**
**  4) On successful return, the array pv contains the following:
**
**        pv[0][0]   x      }
**        pv[0][1]   y      } heliocentric position, au
**        pv[0][2]   z      }
**
**        pv[1][0]   xdot   }
**        pv[1][1]   ydot   } heliocentric velocity, au/d
**        pv[1][2]   zdot   }
**
**     The reference frame is equatorial and is with respect to the
**     mean equator and equinox of epoch J2000.0.
**
**  5) The algorithm is due to J.L. Simon, P. Bretagnon, J. Chapront,
**     M. Chapront-Touze, G. Francou and J. Laskar (Bureau des
**     Longitudes, Paris, France).  From comparisons with JPL
**     ephemeris DE102, they quote the following maximum errors
**     over the interval 1800-2050:
**
**                     L (arcsec)    B (arcsec)      R (km)
**
**        Mercury          4             1             300
**        Venus            5             1             800
**        EMB              6             1            1000
**        Mars            17             1            7700
**        Jupiter         71             5           76000
**        Saturn          81            13          267000
**        Uranus          86             7          712000
**        Neptune         11             1          253000
**
**     Over the interval 1000-3000, they report that the accuracy is no
**     worse than 1.5 times that over 1800-2050.  Outside 1000-3000 the
**     accuracy declines.
**
**     Comparisons of the present function with the JPL DE200 ephemeris
**     give the following RMS errors over the interval 1960-2025:
**
**                      position (km)     velocity (m/s)
**
**        Mercury            334               0.437
**        Venus             1060               0.855
**        EMB               2010               0.815
**        Mars              7690               1.98
**        Jupiter          71700               7.70
**        Saturn          199000              19.4
**        Uranus          564000              16.4
**        Neptune         158000              14.4
**
**     Comparisons against DE200 over the interval 1800-2100 gave the
**     following maximum absolute differences.  (The results using
**     DE406 were essentially the same.)
**
**                   L (arcsec)   B (arcsec)     R (km)   Rdot (m/s)
**
**        Mercury        7            1            500       0.7
**        Venus          7            1           1100       0.9
**        EMB            9            1           1300       1.0
**        Mars          26            1           9000       2.5
**        Jupiter       78            6          82000       8.2
**        Saturn        87           14         263000      24.6
**        Uranus        86            7         661000      27.4
**        Neptune       11            2         248000      21.4
**
**  6) The present SOFA re-implementation of the original Simon et al.
**     Fortran code differs from the original in the following respects:
**
**       *  C instead of Fortran.
**
**       *  The date is supplied in two parts.
**
**       *  The result is returned only in equatorial Cartesian form;
**          the ecliptic longitude, latitude and radius vector are not
**          returned.
**
**       *  The result is in the J2000.0 equatorial frame, not ecliptic.
**
**       *  More is done in-line: there are fewer calls to subroutines.
**
**       *  Different error/warning status values are used.
**
**       *  A different Kepler's-equation-solver is used (avoiding
**          use of double precision complex).
**
**       *  Polynomials in t are nested to minimize rounding errors.
**
**       *  Explicit double constants are used to avoid mixed-mode
**          expressions.
**
**     None of the above changes affects the result significantly.
**
**  7) The returned status indicates the most serious condition
**     encountered during execution of the function.  Illegal np is
**     considered the most serious, overriding failure to converge,
**     which in turn takes precedence over the remote date warning.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**
**  Reference:  Simon, J.L, Bretagnon, P., Chapront, J.,
**              Chapront-Touze, M., Francou, G., and Laskar, J.,
**              Astron.Astrophys., 282, 663 (1994).
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    pv = zeros(shape = (2,3), dtype  = float)
    s = _sofa.iauPlan94(date1,date2,np,pv)
    if s != 0:
        warnings.warn(plan94_msg[s], UserWarning, 2)
    return pv


################################################################
#####    Precession, nutation, polar motion
###############################################################
#iauBi00
try:
    _sofa.iauBi00.argtypes = [
        POINTER(c_double),#dpsibi
        POINTER(c_double),#depsbi
        POINTER(c_double)]#dra
except AttributeError:
    pass
def bi00():
    """
    Frame bias components of IAU 2000 precession-nutation models (part of MHB2000 with additions).

    Returns
    -------
    dpsibi
        longitude and obliquity corrections.
    depsbi
        longitude and obliquity corrections
    dra
        the ICRS RA of the J2000.0 mean equinox.
/*
**  - - - - - - - -
**   i a u B i 0 0
**  - - - - - - - -
**
**  Frame bias components of IAU 2000 precession-nutation models (part
**  of MHB2000 with additions).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Returned:
**     dpsibi,depsbi  double  longitude and obliquity corrections
**     dra            double  the ICRS RA of the J2000.0 mean equinox
**
**  Notes:
**
**  1) The frame bias corrections in longitude and obliquity (radians)
**     are required in order to correct for the offset between the GCRS
**     pole and the mean J2000.0 pole.  They define, with respect to the
**     GCRS frame, a J2000.0 mean pole that is consistent with the rest
**     of the IAU 2000A precession-nutation model.
**
**  2) In addition to the displacement of the pole, the complete
**     description of the frame bias requires also an offset in right
**     ascension.  This is not part of the IAU 2000A model, and is from
**     Chapront et al. (2002).  It is returned in radians.
**
**  3) This is a supplemented implementation of one aspect of the IAU
**     2000A nutation model, formally adopted by the IAU General
**     Assembly in 2000, namely MHB2000 (Mathews et al. 2002).
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G., Astron.
**     Astrophys., 387, 700, 2002.
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
**     and precession   New nutation series for nonrigid Earth and
**     insights into the Earth's interior", J.Geophys.Res., 107, B4,
**     2002.  The MHB2000 code itself was obtained on 9th September 2002
**     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsibi = c_double()
    depsbi = c_double()
    dra = c_double()
    _sofa.iauBi00(byref(dpsibi), byref(depsbi), byref(dra))
    return dpsibi.value, depsbi.value, dra.value


#iauBp00
try:
    _sofa.iauBp00.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float),#rb
        ndpointer(shape=(3,3),dtype = float),#rp
        ndpointer(shape=(3,3),dtype = float)]#rbp
except  AttributeError:
    pass
def bp00(date1,date2):
    """
Frame bias and precession, IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
       TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rb : double[3][3]
        frame bias matrix (Note 2).
    rp : double[3][3]
        precession matrix (Note 3).
    rbp : double[3][3]
        bias-precession matrix (Note 4).
    /*
**  - - - - - - - -
**   i a u B p 0 0
**  - - - - - - - -
**
**  Frame bias and precession, IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rb           double[3][3]   frame bias matrix (Note 2)
**     rp           double[3][3]   precession matrix (Note 3)
**     rbp          double[3][3]   bias-precession matrix (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             date1         date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
**     applying frame bias.
**
**  3) The matrix rp transforms vectors from J2000.0 mean equator and
**     equinox to mean equator and equinox of date by applying
**     precession.
**
**  4) The matrix rbp transforms vectors from GCRS to mean equator and
**     equinox of date by applying frame bias then precession.  It is
**     the product rp x rb.
**
**  5) It is permissible to re-use the same array in the returned
**     arguments.  The arrays are filled in the order given.
**
**  Called:
**     iauBi00      frame bias components, IAU 2000
**     iauPr00      IAU 2000 precession adjustments
**     iauIr        initialize r-matrix to identity
**     iauRx        rotate around X-axis
**     iauRy        rotate around Y-axis
**     iauRz        rotate around Z-axis
**     iauCr        copy r-matrix
**     iauRxr       product of two r-matrices
**
**  Reference:
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  This revision:  2013 August 21
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/

    """
    rb = zeros(shape = (3,3), dtype = float)
    rp = zeros(shape = (3,3), dtype = float)
    rbp = zeros(shape = (3,3), dtype = float)
    _sofa.iauBp00(date1,date2,rb,rp,rbp)
    return rb,rp,rbp


#iauBp06
try:
    _sofa.iauBp06.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float),#rb
        ndpointer(shape=(3,3),dtype = float),#rp
        ndpointer(shape=(3,3),dtype = float)]#rbp
except  AttributeError:
    pass
def bp06(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rb : double[3][3]
        frame bias matrix (Note 2).
    rp : double[3][3]
        precession matrix (Note 3).
    rbp : double[3][3]
        bias-precession matrix (Note 4).
/*
**  - - - - - - - -
**   i a u B p 0 6
**  - - - - - - - -
**
**  Frame bias and precession, IAU 2006.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rb           double[3][3]   frame bias matrix (Note 2)
**     rp           double[3][3]   precession matrix (Note 3)
**     rbp          double[3][3]   bias-precession matrix (Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**             date1         date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rb transforms vectors from GCRS to mean J2000.0 by
**     applying frame bias.
**
**  3) The matrix rp transforms vectors from mean J2000.0 to mean of
**     date by applying precession.
**
**  4) The matrix rbp transforms vectors from GCRS to mean of date by
**     applying frame bias then precession.  It is the product rp x rb.
**
**  5) It is permissible to re-use the same array in the returned
**     arguments.  The arrays are filled in the order given.
**
**  Called:
**     iauPfw06     bias-precession F-W angles, IAU 2006
**     iauFw2m      F-W angles to r-matrix
**     iauPmat06    PB matrix, IAU 2006
**     iauTr        transpose r-matrix
**     iauRxr       product of two r-matrices
**     iauCr        copy r-matrix
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 August 21
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rb = zeros(shape = (3,3), dtype = float)
    rp = zeros(shape = (3,3), dtype = float)
    rbp = zeros(shape = (3,3), dtype = float)
    _sofa.iauBp06(date1,date2,rb,rp,rbp)
    return rb,rp,rbp


#iauBpn2xy
try:
    _sofa.iauBpn2xy.argtypes = [
        ndpointer(shape=(3,3),dtype = float),#rbpn
        POINTER(c_double),#x
        POINTER(c_double)]#y
except  AttributeError:
    pass
def bpn2xy(rbpn):
    """
    Extract from the bias-precession-nutation matrix the X,Y coordinates of the Celestial Intermediate Pole.

    Parameters
    ----------
    rbpn : double[3][3]
        celestial-to-true matrix (Note 1).

    Returns
    -------
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).
/*
**  - - - - - - - - - -
**   i a u B p n 2 x y
**  - - - - - - - - - -
**
**  Extract from the bias-precession-nutation matrix the X,Y coordinates
**  of the Celestial Intermediate Pole.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rbpn      double[3][3]  celestial-to-true matrix (Note 1)
**
**  Returned:
**     x,y       double        Celestial Intermediate Pole (Note 2)
**
**  Notes:
**
**  1) The matrix rbpn transforms vectors from GCRS to true equator (and
**     CIO or equinox) of date, and therefore the Celestial Intermediate
**     Pole unit vector is the bottom row of the matrix.
**
**  2) The arguments x,y are components of the Celestial Intermediate
**     Pole unit vector in the Geocentric Celestial Reference System.
**
**  Reference:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    _sofa.iauBpn2xy(rbpn,byref(x),byref(y))
    return x.value, y.value


#iauC2i00a
try:
    _sofa.iauC2i00a.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2i00a(date1,date2):
    """
    Form the celestial-to-intermediate matrix for a given date using the IAU 2000A precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rc2i : double[3][3]
        celestial-to-intermediate matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u C 2 i 0 0 a
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2000A precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the iauC2i00b function.
**
**  Called:
**     iauPnm00a    classical NPB matrix, IAU 2000A
**     iauC2ibpn    celestial-to-intermediate matrix, given NPB matrix
**
**  References:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2i00a(date1,date2,rc2i)
    return rc2i
    

#iauC2i00b
try:
    _sofa.iauC2i00b.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2i00b(date1,date2):
    """
    Form the celestial-to-intermediate matrix for a given date using the IAU 2000B precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rc2i : double[3][3]
        celestial-to-intermediate matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u C 2 i 0 0 b
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2000B precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the iauC2i00a function.
**
**  Called:
**     iauPnm00b    classical NPB matrix, IAU 2000B
**     iauC2ibpn    celestial-to-intermediate matrix, given NPB matrix
**
**  References:
**
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154
**     (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2i00b(date1,date2,rc2i)
    return rc2i


#iauC2i06a
try:
    _sofa.iauC2i06a.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2i06a(date1,date2):
    """
    Form the celestial-to-intermediate matrix for a given date using the IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rc2i : double[3][3]
        celestial-to-intermediate matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u C 2 i 0 6 a
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date using the
**  IAU 2006 precession and IAU 2000A nutation models.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS]  =  RPOM * R_3(ERA) * rc2i * [CRS]
**
**               =  RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  Called:
**     iauPnm06a    classical NPB matrix, IAU 2006/2000A
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS06       the CIO locator s, given X,Y, IAU 2006
**     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2i06a(date1,date2,rc2i)
    return rc2i

#iauC2ibpn
try:
    _sofa.iauC2ibpn.argtypes = [
        c_double,#date1
        c_double,#date2
        ndpointer(shape=(3,3),dtype = float),#rbpn
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2ibpn(date1,date2,rbpn):
    """
    Form the celestial-to-intermediate matrix for a given date given the bias-precession-nutation matrix.  IAU 2000.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    rbpn : double[3][3]
        celestial-to-true matrix (Note 2).

    Returns
    -------
    rc2i : double[3][3] 
        celestial-to-intermediate matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u C 2 i b p n
**  - - - - - - - - - -
**
**  Form the celestial-to-intermediate matrix for a given date given
**  the bias-precession-nutation matrix.  IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**     rbpn        double[3][3] celestial-to-true matrix (Note 2)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix rbpn transforms vectors from GCRS to true equator (and
**     CIO or equinox) of date.  Only the CIP (bottom row) is used.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  4) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauC2ixy     celestial-to-intermediate matrix, given X,Y
**
**  References:
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2ibpn(date1,date2,rbpn,rc2i)
    return rc2i


#iauC2ixy
try:
    _sofa.iauC2ixy.argtypes = [
        c_double,#date1
        c_double,#date2
        c_double,# x
        c_double,#y
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2ixy(date1,date2,x,y):
    """
Form the celestial to intermediate-frame-of-date matrix for a given date when the CIP X,Y coordinates are known.  IAU 2000.    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).

    Returns
    -------
    rc2i : double[3][3]
        celestial-to-intermediate matrix (Note 3).
/*
**  - - - - - - - - -
**   i a u C 2 i x y
**  - - - - - - - - -
**
**  Form the celestial to intermediate-frame-of-date matrix for a given
**  date when the CIP X,Y coordinates are known.  IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**     x,y         double       Celestial Intermediate Pole (Note 2)
**
**  Returned:
**     rc2i        double[3][3] celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y components
**     of the unit vector in the Geocentric Celestial Reference System.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  4) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     iauC2ixys    celestial-to-intermediate matrix, given X,Y and s
**     iauS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2ixy(date1,date2,x,y,rc2i)
    return rc2i


#iauC2ixys
try:
    _sofa.iauC2ixys.argtypes = [
        c_double,# x
        c_double,#y
        c_double,#s
        ndpointer(shape=(3,3),dtype = float)]#rc2i
except  AttributeError:
    pass
def c2ixys(x,y,s):
    """
  Form the celestial to intermediate-frame-of-date matrix given the CIP X,Y and the CIO locator s.  

    Parameters
    ----------
    x : float
        Celestial Intermediate Pole (Note 1).
    y : float
        Celestial Intermediate Pole (Note 1).
    s : float
        the CIO locator s (Note 2).

    Returns
    -------
    rc2i : float[3][3]
        celestial-to-intermediate matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u C 2 i x y s
**  - - - - - - - - - -
**
**  Form the celestial to intermediate-frame-of-date matrix given the CIP
**  X,Y and the CIO locator s.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     x,y      double         Celestial Intermediate Pole (Note 1)
**     s        double         the CIO locator s (Note 2)
**
**  Returned:
**     rc2i     double[3][3]   celestial-to-intermediate matrix (Note 3)
**
**  Notes:
**
**  1) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  2) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  3) The matrix rc2i is the first stage in the transformation from
**     celestial to terrestrial coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = RC2T * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  Called:
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRy        rotate around Y-axis
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2014 November 7
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2i = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2ixys(x,y,s,rc2i)
    return rc2i


#iauC2t00a
try:
    _sofa.iauC2t00a.argtypes = [
        c_double,# tta
        c_double,# ttb
        c_double,# uta
        c_double,# utb
        c_double,# xp
        c_double,# yp
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2t00a(tta,ttb,uta,utb,xp,yp):
    """
  Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000A nutation model.  

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date (Note 1).
    ttb : float
        TT as a 2-part Julian Date (Note 1).
    uta : float
        UT1 as a 2-part Julian Date (Note 1).
    utb : float
        UT1 as a 2-part Julian Date (Note 1).
    xp : float
        coordinates of the pole (radians, Note 2).
    yp : float
        coordinates of the pole (radians, Note 2).

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u C 2 t 0 0 a
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2000A nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  4) A faster, but slightly less accurate result (about 1 mas), can
**     be obtained by using instead the iauC2t00b function.
**
**  Called:
**     iauC2i00a    celestial-to-intermediate matrix, IAU 2000A
**     iauEra00     Earth rotation angle, IAU 2000
**     iauSp00      the TIO locator s', IERS 2000
**     iauPom00     polar motion matrix
**     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2t00a(tta,ttb,uta,utb,xp,yp,rc2t)
    return rc2t


#iauC2t00b
try:
    _sofa.iauC2t00b.argtypes = [
        c_double,# tta
        c_double,# ttb
        c_double,# uta
        c_double,# utb
        c_double,# xp
        c_double,# yp
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2t00b(tta,ttb,uta,utb,xp,yp):
    """
    Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2000B nutation model.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date (Note 1).
    ttb : float
        TT as a 2-part Julian Date (Note 1).
    uta : float
        UT1 as a 2-part Julian Date (Note 1).
    utb : float
        UT1 as a 2-part Julian Date (Note 1).
    xp : float
        coordinates of the pole (radians, Note 2).
    yp : float
        coordinates of the pole (radians, Note 2).

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u C 2 t 0 0 b
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2000B nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  4) The present function is faster, but slightly less accurate (about
**     1 mas), than the iauC2t00a function.
**
**  Called:
**     iauC2i00b    celestial-to-intermediate matrix, IAU 2000B
**     iauEra00     Earth rotation angle, IAU 2000
**     iauPom00     polar motion matrix
**     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2t00b(tta,ttb,uta,utb,xp,yp,rc2t)
    return rc2t


#iauC2t06a
try:
    _sofa.iauC2t06a.argtypes = [
        c_double,# tta
        c_double,# ttb
        c_double,# uta
        c_double,# utb
        c_double,# xp
        c_double,# yp
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2t06a(tta,ttb,uta,utb,xp,yp):
    """
    Form the celestial to terrestrial matrix given the date, the UT1 and the polar motion, using the IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date (Note 1).
    ttb : float
        TT as a 2-part Julian Date (Note 1).
    uta : float
        UT1 as a 2-part Julian Date (Note 1).
    utb : float
        UT1 as a 2-part Julian Date (Note 1).
    xp : float
        coordinates of the pole (radians, Note 2).
    yp : float
        coordinates of the pole (radians, Note 2).

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u C 2 t 0 6 a
**  - - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1 and
**  the polar motion, using the IAU 2006 precession and IAU 2000A
**  nutation models.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     xp,yp    double         coordinates of the pole (radians, Note 2)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 3)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  3) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RC2I is the
**     celestial-to-intermediate matrix, ERA is the Earth rotation
**     angle and RPOM is the polar motion matrix.
**
**  Called:
**     iauC2i06a    celestial-to-intermediate matrix, IAU 2006/2000A
**     iauEra00     Earth rotation angle, IAU 2000
**     iauSp00      the TIO locator s', IERS 2000
**     iauPom00     polar motion matrix
**     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2t06a(tta,ttb,uta,utb,xp,yp,rc2t)
    return rc2t


#iauC2tcio
try:
    _sofa.iauC2tcio.argtypes = [
        ndpointer(shape=(3,3),dtype = float),#rc2i
        c_double,# era
        ndpointer(shape=(3,3),dtype = float),#rpom
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2tcio(rc2i,era,rpom):
    """
    Assemble the celestial to terrestrial matrix from CIO-based components (the celestial-to-intermediate matrix, the Earth Rotation  Angle and the polar motion matrix).

    Parameters
    ----------
    rc2i : double[3][3]
        celestial-to-intermediate matrix.
    era : float
        Earth rotation angle (radians).
    rpom : double[3][3]
        polar-motion matrix.

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix.
/*
**  - - - - - - - - - -
**   i a u C 2 t c i o
**  - - - - - - - - - -
**
**  Assemble the celestial to terrestrial matrix from CIO-based
**  components (the celestial-to-intermediate matrix, the Earth Rotation
**  Angle and the polar motion matrix).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rc2i     double[3][3]    celestial-to-intermediate matrix
**     era      double          Earth rotation angle (radians)
**     rpom     double[3][3]    polar-motion matrix
**
**  Returned:
**     rc2t     double[3][3]    celestial-to-terrestrial matrix
**
**  Notes:
**
**  1) This function constructs the rotation matrix that transforms
**     vectors in the celestial system into vectors in the terrestrial
**     system.  It does so starting from precomputed components, namely
**     the matrix which rotates from celestial coordinates to the
**     intermediate frame, the Earth rotation angle and the polar motion
**     matrix.  One use of the present function is when generating a
**     series of celestial-to-terrestrial matrices where only the Earth
**     Rotation Angle changes, avoiding the considerable overhead of
**     recomputing the precession-nutation more often than necessary to
**     achieve given accuracy objectives.
**
**  2) The relationship between the arguments is as follows:
**
**        [TRS] = RPOM * R_3(ERA) * rc2i * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003).
**
**  Called:
**     iauCr        copy r-matrix
**     iauRz        rotate around Z-axis
**     iauRxr       product of two r-matrices
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2013 August 24
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2tcio(rc2i,era,rpom,rc2t)
    return rc2t


#iauC2teqx
try:
    _sofa.iauC2teqx.argtypes = [
        ndpointer(shape=(3,3),dtype = float),#rbpn
        c_double,# gst
        ndpointer(shape=(3,3),dtype = float),#rpom
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2teqx(rbpn,gst,rpom):
    """
    

    Parameters
    ----------
    rbpn : double[3][3]
        celestial-to-true matrix.
    gst : float
        Greenwich (apparent) Sidereal Time (radians).
    rpom : double[3][3]
        polar-motion matrix.

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u C 2 t e q x
**  - - - - - - - - - -
**
**  Assemble the celestial to terrestrial matrix from equinox-based
**  components (the celestial-to-true matrix, the Greenwich Apparent
**  Sidereal Time and the polar motion matrix).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rbpn   double[3][3]  celestial-to-true matrix
**     gst    double        Greenwich (apparent) Sidereal Time (radians)
**     rpom   double[3][3]  polar-motion matrix
**
**  Returned:
**     rc2t   double[3][3]  celestial-to-terrestrial matrix (Note 2)
**
**  Notes:
**
**  1) This function constructs the rotation matrix that transforms
**     vectors in the celestial system into vectors in the terrestrial
**     system.  It does so starting from precomputed components, namely
**     the matrix which rotates from celestial coordinates to the
**     true equator and equinox of date, the Greenwich Apparent Sidereal
**     Time and the polar motion matrix.  One use of the present function
**     is when generating a series of celestial-to-terrestrial matrices
**     where only the Sidereal Time changes, avoiding the considerable
**     overhead of recomputing the precession-nutation more often than
**     necessary to achieve given accuracy objectives.
**
**  2) The relationship between the arguments is as follows:
**
**        [TRS] = rpom * R_3(gst) * rbpn * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003).
**
**  Called:
**     iauCr        copy r-matrix
**     iauRz        rotate around Z-axis
**     iauRxr       product of two r-matrices
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 August 24
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2teqx(rbpn,gst,rpom,rc2t)
    return rc2t


#iauC2tpe
try:
    _sofa.iauC2tpe.argtypes = [
        c_double,# tta
        c_double,# ttb
        c_double,# uta
        c_double,# utb
        c_double,# dpsi
        c_double,# deps
        c_double,# xp
        c_double,# yp
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp):
    """
    

    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date (Note 1).
    ttb : float
        TT as a 2-part Julian Date (Note 1).
    uta : float
        UT1 as a 2-part Julian Date (Note 1).
    utb : float
        UT1 as a 2-part Julian Date (Note 1).
    dpsi : float
        nutation (Note 2).
    deps : float
        nutation (Note 2).
    xp : float
        coordinates of the pole (radians, Note 3).
    yp : float
        coordinates of the pole (radians, Note 3).

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 4).
/*
**  - - - - - - - - -
**   i a u C 2 t p e
**  - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1,
**  the nutation and the polar motion.  IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     tta,ttb    double        TT as a 2-part Julian Date (Note 1)
**     uta,utb    double        UT1 as a 2-part Julian Date (Note 1)
**     dpsi,deps  double        nutation (Note 2)
**     xp,yp      double        coordinates of the pole (radians, Note 3)
**
**  Returned:
**     rc2t       double[3][3]  celestial-to-terrestrial matrix (Note 4)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any of
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.  For high-accuracy
**     applications, free core nutation should be included as well as
**     any other relevant corrections to the position of the CIP.
**
**  3) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  4) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(GST) * RBPN * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), RBPN is the
**     bias-precession-nutation matrix, GST is the Greenwich (apparent)
**     Sidereal Time and RPOM is the polar motion matrix.
**
**  5) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     iauPn00      bias/precession/nutation results, IAU 2000
**     iauGmst00    Greenwich mean sidereal time, IAU 2000
**     iauSp00      the TIO locator s', IERS 2000
**     iauEe00      equation of the equinoxes, IAU 2000
**     iauPom00     polar motion matrix
**     iauC2teqx    form equinox-based celestial-to-terrestrial matrix
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2tpe(tta,ttb,uta,utb,dpsi,deps,xp,yp,rc2t)
    return rc2t


#iauC2txy
try:
    _sofa.iauC2txy.argtypes = [
        c_double,# tta
        c_double,# ttb
        c_double,# uta
        c_double,# utb
        c_double,# x
        c_double,# y
        c_double,# xp
        c_double,# yp
        ndpointer(shape=(3,3),dtype = float)]#rc2t
except  AttributeError:
    pass
def c2txy(tta,ttb,uta,utb,x,y,xp,yp):
    """
    
Form the celestial to terrestrial matrix given the date, the UT1, the CIP coordinates and the polar motion.  IAU 2000.
    Parameters
    ----------
    tta : float
        TT as a 2-part Julian Date (Note 1).
    ttb : float
        TT as a 2-part Julian Date (Note 1).
    uta : float
        UT1 as a 2-part Julian Date (Note 1).
    utb : float
        UT1 as a 2-part Julian Date (Note 1).
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).
    xp : float
        coordinates of the pole (radians, Note 3).
    yp : float
        coordinates of the pole (radians, Note 3).

    Returns
    -------
    rc2t : double[3][3]
        celestial-to-terrestrial matrix (Note 4).
/*
**  - - - - - - - - -
**   i a u C 2 t x y
**  - - - - - - - - -
**
**  Form the celestial to terrestrial matrix given the date, the UT1,
**  the CIP coordinates and the polar motion.  IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     tta,ttb  double         TT as a 2-part Julian Date (Note 1)
**     uta,utb  double         UT1 as a 2-part Julian Date (Note 1)
**     x,y      double         Celestial Intermediate Pole (Note 2)
**     xp,yp    double         coordinates of the pole (radians, Note 3)
**
**  Returned:
**     rc2t     double[3][3]   celestial-to-terrestrial matrix (Note 4)
**
**  Notes:
**
**  1) The TT and UT1 dates tta+ttb and uta+utb are Julian Dates,
**     apportioned in any convenient way between the arguments uta and
**     utb.  For example, JD(UT1)=2450123.7 could be expressed in any o
**     these ways, among others:
**
**             uta            utb
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 and MJD methods are good compromises
**     between resolution and convenience.  In the case of uta,utb, the
**     date & time method is best matched to the Earth rotation angle
**     algorithm used:  maximum precision is delivered when the uta
**     argument is for 0hrs UT1 on the day in question and the utb
**     argument lies in the range 0 to 1, or vice versa.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  4) The matrix rc2t transforms from celestial to terrestrial
**     coordinates:
**
**        [TRS] = RPOM * R_3(ERA) * RC2I * [CRS]
**
**              = rc2t * [CRS]
**
**     where [CRS] is a vector in the Geocentric Celestial Reference
**     System and [TRS] is a vector in the International Terrestrial
**     Reference System (see IERS Conventions 2003), ERA is the Earth
**     Rotation Angle and RPOM is the polar motion matrix.
**
**  5) Although its name does not include "00", This function is in fact
**     specific to the IAU 2000 models.
**
**  Called:
**     iauC2ixy     celestial-to-intermediate matrix, given X,Y
**     iauEra00     Earth rotation angle, IAU 2000
**     iauSp00      the TIO locator s', IERS 2000
**     iauPom00     polar motion matrix
**     iauC2tcio    form CIO-based celestial-to-terrestrial matrix
**
** Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rc2t = zeros(shape=(3,3),dtype = float)
    _sofa.iauC2txy(tta,ttb,uta,utb,x,y,xp,yp,rc2t)
    return rc2t


#iauEo06
try:
    _sofa.iauEo06a.argtypes = [
        c_double,# date1
        c_double,# date2
        ]
    _sofa.iauEo06a.restype = c_double #eo
except  AttributeError:
    pass
def eo06a(date1, date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    eo : float
        equation of the origins in radians.
/*
**  - - - - - - - - -
**   i a u E o 0 6 a
**  - - - - - - - - -
**
**  Equation of the origins, IAU 2006 precession and IAU 2000A nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    equation of the origins in radians
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The equation of the origins is the distance between the true
**     equinox and the celestial intermediate origin and, equivalently,
**     the difference between Earth rotation angle and Greenwich
**     apparent sidereal time (ERA-GST).  It comprises the precession
**     (since J2000.0) in right ascension plus the equation of the
**     equinoxes (including the small correction terms).
**
**  Called:
**     iauPnm06a    classical NPB matrix, IAU 2006/2000A
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS06       the CIO locator s, given X,Y, IAU 2006
**     iauEors      equation of the origins, given NPB matrix and s
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eo = _sofa.iauEo06a(date1, date2)
    return eo

#iauEors
try:
    _sofa.iauEors.argtypes = [
        ndpointer(shape = (3,3),dtype = float),# rnpb[3][3]
        c_double,# s
        ]
    _sofa.iauEors.restype = c_double #eo
except  AttributeError:
    pass
def eors(rnpb,s):
    """
    Equation of the origins, given the classical NPB matrix and the quantity s.

    Parameters
    ----------
    rnpb : double[3][3] 
        classical nutation x precession x bias matrix.
    s : float
        the quantity s (the CIO locator).

    Returns
    -------
    eo : float
        the equation of the origins in radians..
/*
**  - - - - - - - -
**   i a u E o r s
**  - - - - - - - -
**
**  Equation of the origins, given the classical NPB matrix and the
**  quantity s.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rnpb  double[3][3]  classical nutation x precession x bias matrix
**     s     double        the quantity s (the CIO locator)
**
**  Returned (function value):
**           double        the equation of the origins in radians.
**
**  Notes:
**
**  1)  The equation of the origins is the distance between the true
**      equinox and the celestial intermediate origin and, equivalently,
**      the difference between Earth rotation angle and Greenwich
**      apparent sidereal time (ERA-GST).  It comprises the precession
**      (since J2000.0) in right ascension plus the equation of the
**      equinoxes (including the small correction terms).
**
**  2)  The algorithm is from Wallace & Capitaine (2006).
**
** References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eo = _sofa.iauEors(rnpb,s)
    return eo


#iauFw2m
try:
    _sofa.iauFw2m.argtypes = [
        c_double, #gamb
        c_double, #phib
        c_double, #psi,
        c_double, #eps
        ndpointer(shape = (3,3), dtype =float),# r[3][3]
        ]
except  AttributeError:
    pass
def fw2m(gamb,phib,psi,eps):
    """
    Form rotation matrix given the Fukushima-Williams angles.

    Parameters
    ----------
    gamb : float
        F-W angle gamma_bar (radians).
    phib : float
        F-W angle phi_bar (radians).
    psi : float
        F-W angle psi (radians).
    eps : float
        F-W angle epsilon (radians).

    Returns
    -------
    r : double[3][3]
        rotation matrix.
/*
**  - - - - - - - -
**   i a u F w 2 m
**  - - - - - - - -
**
**  Form rotation matrix given the Fukushima-Williams angles.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     gamb     double         F-W angle gamma_bar (radians)
**     phib     double         F-W angle phi_bar (radians)
**     psi      double         F-W angle psi (radians)
**     eps      double         F-W angle epsilon (radians)
**
**  Returned:
**     r        double[3][3]   rotation matrix
**
**  Notes:
**
**  1) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = ecliptic pole of date,
**     and   P = CIP,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma = epE
**        phib = phi = pE
**        psi = psi = pEP
**        eps = epsilon = EP
**
**  2) The matrix representing the combined effects of frame bias,
**     precession and nutation is:
**
**        NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)
**
**  3) Three different matrices can be constructed, depending on the
**     supplied angles:
**
**     o  To obtain the nutation x precession x frame bias matrix,
**        generate the four precession angles, generate the nutation
**        components and add them to the psi_bar and epsilon_A angles,
**        and call the present function.
**
**     o  To obtain the precession x frame bias matrix, generate the
**        four precession angles and call the present function.
**
**     o  To obtain the frame bias matrix, generate the four precession
**        angles for date J2000.0 and call the present function.
**
**     The nutation-only and precession-only matrices can if necessary
**     be obtained by combining these three appropriately.
**
**  Called:
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRx        rotate around X-axis
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = zeros(shape = (3,3), dtype = float)
    _sofa.iauFw2m(gamb,phib,psi,eps,r)
    return r


#iauFw2xy
try:
    _sofa.iauFw2xy.argtypes = [
        c_double, #gamb
        c_double, #phib
        c_double, #psi,
        c_double, #eps
        POINTER(c_double),# x
        POINTER(c_double),# y
        ]
except  AttributeError:
    pass
def fw2xy(gamb,phib,psi,eps):
    """
    CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

    Parameters
    ----------
    gamb : float
         F-W angle gamma_bar (radians).
    phib : float
        F-W angle phi_bar (radians).
    psi : float
        F-W angle psi (radians).
    eps : float
        F-W angle epsilon (radians).

    Returns
    -------
    x : float
        CIP unit vector X,Y.
    y : float
        CIP unit vector X,Y.
/*
**  - - - - - - - - -
**   i a u F w 2 x y
**  - - - - - - - - -
**
**  CIP X,Y given Fukushima-Williams bias-precession-nutation angles.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     gamb     double    F-W angle gamma_bar (radians)
**     phib     double    F-W angle phi_bar (radians)
**     psi      double    F-W angle psi (radians)
**     eps      double    F-W angle epsilon (radians)
**
**  Returned:
**     x,y      double    CIP unit vector X,Y
**
**  Notes:
**
**  1) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole
**           E = ecliptic pole of date,
**     and   P = CIP,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma = epE
**        phib = phi = pE
**        psi = psi = pEP
**        eps = epsilon = EP
**
**  2) The matrix representing the combined effects of frame bias,
**     precession and nutation is:
**
**        NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)
**
**     The returned values x,y are elements [2][0] and [2][1] of the
**     matrix.  Near J2000.0, they are essentially angles in radians.
**
**  Called:
**     iauFw2m      F-W angles to r-matrix
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  This revision:  2013 September 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    _sofa.iauFw2xy(gamb,phib,psi,eps,byref(x),byref(y))
    return x.value, y.value


#iauLtp
try:
    _sofa.iauLtp.argtypes = [
        c_double, #epj
        ndpointer(shape = (3,3), dtype = float)#rp[3][3]
        ]
except  AttributeError:
    pass
def ltp(epj):
    """
    Long-term precession matrix.

    Parameters
    ----------
    epj : float
        Julian epoch (TT).

    Returns
    -------
    rp : double[3][3]
        precession matrix, J2000.0 to date.
/*
**  - - - - - - -
**   i a u L t p
**  - - - - - - -
**
**  Long-term precession matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     rp      double[3][3]   precession matrix, J2000.0 to date
**
**  Notes:
**
**  1) The matrix is in the sense
**
**        P_date = rp x P_J2000,
**
**     where P_J2000 is a vector with respect to the J2000.0 mean
**     equator and equinox and P_date is the same vector with respect to
**     the equator and equinox of epoch epj.
**
**  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  Called:
**     iauLtpequ    equator pole, long term
**     iauLtpecl    ecliptic pole, long term
**     iauPxp       vector product
**     iauPn        normalize vector
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2015 December 6
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rp = zeros(shape = (3,3), dtype = float)
    _sofa.iauLtp(epj,rp)
    return rp


#iauLtpb
try:
    _sofa.iauLtpb.argtypes = [
        c_double, #epj
        ndpointer(shape = (3,3), dtype = float)#rpb[3][3]
        ]
except  AttributeError:
    pass
def ltpb(epj):
    """
    Long-term precession matrix, including ICRS frame bias.

    Parameters
    ----------
    epj : float
         Julian epoch (TT).

    Returns
    -------
    rpb : double[3][3]
        precession-bias matrix, J2000.0 to date.
/*
**  - - - - - - - -
**   i a u L t p b
**  - - - - - - - -
**
**  Long-term precession matrix, including ICRS frame bias.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     rpb     double[3][3]   precession-bias matrix, J2000.0 to date
**
**  Notes:
**
**  1) The matrix is in the sense
**
**        P_date = rpb x P_ICRS,
**
**     where P_ICRS is a vector in the Geocentric Celestial Reference
**     System, and P_date is the vector with respect to the Celestial
**     Intermediate Reference System at that date but with nutation
**     neglected.
**
**  2) A first order frame bias formulation is used, of sub-
**     microarcsecond accuracy compared with a full 3D rotation.
**
**  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2015 December 6
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rpb = zeros(shape = (3,3), dtype = float)
    _sofa.iauLtpb(epj,rpb)
    return rpb


#iauLtpecl
try:
    _sofa.iauLtpecl.argtypes = [
        c_double, #epj
        ndpointer(shape = 3, dtype = float)#vec[3]
        ]
except  AttributeError:
    pass
def ltpecl(epj):
    """
    Long-term precession of the ecliptic.

    Parameters
    ----------
    epj : float
        Julian epoch (TT).

    Returns
    -------
    vec : float[3]
        ecliptic pole unit vector.
/*
**  - - - - - - - - - -
**   i a u L t p e c l
**  - - - - - - - - - -
**
**  Long-term precession of the ecliptic.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     vec     double[3]      ecliptic pole unit vector
**
**  Notes:
**
**  1) The returned vector is with respect to the J2000.0 mean equator
**     and equinox.
**
**  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    vec = zeros(shape = 3, dtype = float)
    _sofa.iauLtpecl(epj,vec)
    return vec

#iauLtpequ
try:
    _sofa.iauLtpequ.argtypes = [
        c_double, #epj
        ndpointer(shape = 3, dtype = float)#veq[3][3]
        ]
except  AttributeError:
    pass
def ltpequ(epj):
    """
    Long-term precession of the equator.

    Parameters
    ----------
    epj : float
        Julian epoch (TT).

    Returns
    -------
    veq : float[3]
        equator pole unit vector.
/*
**  - - - - - - - - - -
**   i a u L t p e q u
**  - - - - - - - - - -
**
**  Long-term precession of the equator.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     veq     double[3]      equator pole unit vector
**
**  Notes:
**
**  1) The returned vector is with respect to the J2000.0 mean equator
**     and equinox.
**
**  2) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    veq = zeros(shape = 3, dtype = float)
    _sofa.iauLtpequ(epj,veq)
    return veq


#iauNum00a
try:
    _sofa.iauNum00a.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape = (3,3), dtype = float) #rmatn[3][3]
        ]
except  AttributeError:
    pass
def num00a(date1,date2):
    """
    Form the matrix of nutation for a given date, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rmatn : double[3][3]
        nutation matrix.
/*
**  - - - - - - - - - -
**   i a u N u m 0 0 a
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2000A model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn        double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the iauNum00b function.
**
**  Called:
**     iauPn00a     bias/precession/nutation, IAU 2000A
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatn = zeros(shape = (3,3), dtype = float)
    _sofa.iauNum00a(date1,date2,rmatn)
    return rmatn 


#iauNum00b
try:
    _sofa.iauNum00b.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape = (3,3), dtype = float) #rmatn[3][3]
        ]
except  AttributeError:
    pass
def num00b(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rmatn : double[3][3]
        nutation matrix.
/*
**  - - - - - - - - - -
**   i a u N u m 0 0 b
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2000B model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn        double[3][3]   nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the iauNum00a function.
**
**  Called:
**     iauPn00b     bias/precession/nutation, IAU 2000B
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatn = zeros(shape = (3,3), dtype = float)
    _sofa.iauNum00b(date1,date2,rmatn)
    return rmatn 


#iauNum06a
try:
    _sofa.iauNum06a.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape = (3,3), dtype = float) #rmatn[3][3]
        ]
except  AttributeError:
    pass
def num06a(date1,date2):
    """
    Form the matrix of nutation for a given date, IAU 2006/2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rmatn : double[3][3]
        nutation matrix.
/*
**  - - - - - - - - - -
**   i a u N u m 0 6 a
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 2006/2000A model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2   double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rmatn         double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(true) = rmatn * V(mean), where
**     the p-vector V(true) is with respect to the true equatorial triad
**     of date and the p-vector V(mean) is with respect to the mean
**     equatorial triad of date.
**
**  Called:
**     iauObl06     mean obliquity, IAU 2006
**     iauNut06a    nutation, IAU 2006/2000A
**     iauNumat     form nutation matrix
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatn = zeros(shape = (3,3), dtype = float)
    _sofa.iauNum06a(date1,date2,rmatn)
    return rmatn 


#iauNumat
try:
    _sofa.iauNumat.argtypes = [
        c_double, #epsa
        c_double, #dpsi
        c_double, #deps
        ndpointer(shape = (3,3), dtype = float) #rmatn[3][3]
        ]
except  AttributeError:
    pass
def numat(epsa,dpsi,deps):
    """
    
Form the matrix of nutation.
    Parameters
    ----------
    epsa : float
        mean obliquity of date (Note 1).
    dpsi : float
        nutation.
    deps : float
        nutation.

    Returns
    -------
    rmatn : double[3][3]
        nutation matrix.
/*
**  - - - - - - - - -
**   i a u N u m a t
**  - - - - - - - - -
**
**  Form the matrix of nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epsa        double         mean obliquity of date (Note 1)
**     dpsi,deps   double         nutation (Note 2)
**
**  Returned:
**     rmatn       double[3][3]   nutation matrix (Note 3)
**
**  Notes:
**
**
**  1) The supplied mean obliquity epsa, must be consistent with the
**     precession-nutation models from which dpsi and deps were obtained.
**
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.
**
**  3) The matrix operates in the sense V(true) = rmatn * V(mean),
**     where the p-vector V(true) is with respect to the true
**     equatorial triad of date and the p-vector V(mean) is with
**     respect to the mean equatorial triad of date.
**
**  Called:
**     iauIr        initialize r-matrix to identity
**     iauRx        rotate around X-axis
**     iauRz        rotate around Z-axis
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222-3 (p114).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatn = zeros(shape = (3,3), dtype = float)
    _sofa.iauNumat(epsa,dpsi,deps,rmatn)
    return rmatn 


#iauNut00a
try:
    _sofa.iauNut00a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        ]
except  AttributeError:
    pass
def nut00a(date1,date2):
    """
    Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation  with free core nutation omitted).

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary (Note 2).
    deps : float
        nutation, luni-solar + planetary (Note 2).
/*
**  - - - - - - - - - -
**   i a u N u t 0 0 a
**  - - - - - - - - - -
**
**  Nutation, IAU 2000A model (MHB2000 luni-solar and planetary nutation
**  with free core nutation omitted).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the equinox and ecliptic of date.  The
**     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
**     value of 84381.448 arcsec.
**
**     Both the luni-solar and planetary nutations are included.  The
**     latter are due to direct planetary nutations and the
**     perturbations of the lunar and terrestrial orbits.
**
**  3) The function computes the MHB2000 nutation series with the
**     associated corrections for planetary nutations.  It is an
**     implementation of the nutation part of the IAU 2000A precession-
**     nutation model, formally adopted by the IAU General Assembly in
**     2000, namely MHB2000 (Mathews et al. 2002), but with the free
**     core nutation (FCN - see Note 4) omitted.
**
**  4) The full MHB2000 model also contains contributions to the
**     nutations in longitude and obliquity due to the free-excitation
**     of the free-core-nutation during the period 1979-2000.  These FCN
**     terms, which are time-dependent and unpredictable, are NOT
**     included in the present function and, if required, must be
**     independently computed.  With the FCN corrections included, the
**     present function delivers a pole which is at current epochs
**     accurate to a few hundred microarcseconds.  The omission of FCN
**     introduces further errors of about that size.
**
**  5) The present function provides classical nutation.  The MHB2000
**     algorithm, from which it is adapted, deals also with (i) the
**     offsets between the GCRS and mean poles and (ii) the adjustments
**     in longitude and obliquity due to the changed precession rates.
**     These additional functions, namely frame bias and precession
**     adjustments, are supported by the SOFA functions iauBi00  and
**     iauPr00.
**
**  6) The MHB2000 algorithm also provides "total" nutations, comprising
**     the arithmetic sum of the frame bias, precession adjustments,
**     luni-solar nutation and planetary nutation.  These total
**     nutations can be used in combination with an existing IAU 1976
**     precession implementation, such as iauPmat76,  to deliver GCRS-
**     to-true predictions of sub-mas accuracy at current dates.
**     However, there are three shortcomings in the MHB2000 model that
**     must be taken into account if more accurate or definitive results
**     are required (see Wallace 2002):
**
**       (i) The MHB2000 total nutations are simply arithmetic sums,
**           yet in reality the various components are successive Euler
**           rotations.  This slight lack of rigor leads to cross terms
**           that exceed 1 mas after a century.  The rigorous procedure
**           is to form the GCRS-to-true rotation matrix by applying the
**           bias, precession and nutation in that order.
**
**      (ii) Although the precession adjustments are stated to be with
**           respect to Lieske et al. (1977), the MHB2000 model does
**           not specify which set of Euler angles are to be used and
**           how the adjustments are to be applied.  The most literal
**           and straightforward procedure is to adopt the 4-rotation
**           epsilon_0, psi_A, omega_A, xi_A option, and to add DPSIPR
**           to psi_A and DEPSPR to both omega_A and eps_A.
**
**     (iii) The MHB2000 model predates the determination by Chapront
**           et al. (2002) of a 14.6 mas displacement between the
**           J2000.0 mean equinox and the origin of the ICRS frame.  It
**           should, however, be noted that neglecting this displacement
**           when calculating star coordinates does not lead to a
**           14.6 mas change in right ascension, only a small second-
**           order distortion in the pattern of the precession-nutation
**           effect.
**
**     For these reasons, the SOFA functions do not generate the "total
**     nutations" directly, though they can of course easily be
**     generated by calling iauBi00, iauPr00 and the present function
**     and adding the results.
**
**  7) The MHB2000 model contains 41 instances where the same frequency
**     appears multiple times, of which 38 are duplicates and three are
**     triplicates.  To keep the present code close to the original MHB
**     algorithm, this small inefficiency has not been corrected.
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFame03    mean longitude of Mercury
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFama03    mean longitude of Mars
**     iauFaju03    mean longitude of Jupiter
**     iauFasa03    mean longitude of Saturn
**     iauFaur03    mean longitude of Uranus
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    _sofa.iauNut00a(date1,date2,byref(dpsi),byref(deps))
    return dpsi.value, deps.value


#iauNut00b
try:
    _sofa.iauNut00b.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        ]
except  AttributeError:
    pass
def nut00b(date1,date2):
    """
    Nutation, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary (Note 2).
    deps : float
        nutation, luni-solar + planetary (Note 2).
/*
**  - - - - - - - - - -
**   i a u N u t 0 0 b
**  - - - - - - - - - -
**
**  Nutation, IAU 2000B model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double    nutation, luni-solar + planetary (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the equinox and ecliptic of date.  The
**     obliquity at J2000.0 is assumed to be the Lieske et al. (1977)
**     value of 84381.448 arcsec.  (The errors that result from using
**     this function with the IAU 2006 value of 84381.406 arcsec can be
**     neglected.)
**
**     The nutation model consists only of luni-solar terms, but
**     includes also a fixed offset which compensates for certain long-
**     period planetary terms (Note 7).
**
**  3) This function is an implementation of the IAU 2000B abridged
**     nutation model formally adopted by the IAU General Assembly in
**     2000.  The function computes the MHB_2000_SHORT luni-solar
**     nutation series (Luzum 2001), but without the associated
**     corrections for the precession rate adjustments and the offset
**     between the GCRS and J2000.0 mean poles.
**
**  4) The full IAU 2000A (MHB2000) nutation model contains nearly 1400
**     terms.  The IAU 2000B model (McCarthy & Luzum 2003) contains only
**     77 terms, plus additional simplifications, yet still delivers
**     results of 1 mas accuracy at present epochs.  This combination of
**     accuracy and size makes the IAU 2000B abridged nutation model
**     suitable for most practical applications.
**
**     The function delivers a pole accurate to 1 mas from 1900 to 2100
**     (usually better than 1 mas, very occasionally just outside
**     1 mas).  The full IAU 2000A model, which is implemented in the
**     function iauNut00a (q.v.), delivers considerably greater accuracy
**     at current dates;  however, to realize this improved accuracy,
**     corrections for the essentially unpredictable free-core-nutation
**     (FCN) must also be included.
**
**  5) The present function provides classical nutation.  The
**     MHB_2000_SHORT algorithm, from which it is adapted, deals also
**     with (i) the offsets between the GCRS and mean poles and (ii) the
**     adjustments in longitude and obliquity due to the changed
**     precession rates.  These additional functions, namely frame bias
**     and precession adjustments, are supported by the SOFA functions
**     iauBi00  and iauPr00.
**
**  6) The MHB_2000_SHORT algorithm also provides "total" nutations,
**     comprising the arithmetic sum of the frame bias, precession
**     adjustments, and nutation (luni-solar + planetary).  These total
**     nutations can be used in combination with an existing IAU 1976
**     precession implementation, such as iauPmat76,  to deliver GCRS-
**     to-true predictions of mas accuracy at current epochs.  However,
**     for symmetry with the iauNut00a  function (q.v. for the reasons),
**     the SOFA functions do not generate the "total nutations"
**     directly.  Should they be required, they could of course easily
**     be generated by calling iauBi00, iauPr00 and the present function
**     and adding the results.
**
**  7) The IAU 2000B model includes "planetary bias" terms that are
**     fixed in size but compensate for long-period nutations.  The
**     amplitudes quoted in McCarthy & Luzum (2003), namely
**     Dpsi = -1.5835 mas and Depsilon = +1.6339 mas, are optimized for
**     the "total nutations" method described in Note 6.  The Luzum
**     (2001) values used in this SOFA implementation, namely -0.135 mas
**     and +0.388 mas, are optimized for the "rigorous" method, where
**     frame bias, precession and nutation are applied separately and in
**     that order.  During the interval 1995-2050, the SOFA
**     implementation delivers a maximum error of 1.001 mas (not
**     including FCN).
**
**  References:
**
**     Lieske, J.H., Lederle, T., Fricke, W., Morando, B., "Expressions
**     for the precession quantities based upon the IAU /1976/ system of
**     astronomical constants", Astron.Astrophys. 58, 1-2, 1-16. (1977)
**
**     Luzum, B., private communication, 2001 (Fortran code
**     MHB_2000_SHORT)
**
**     McCarthy, D.D. & Luzum, B.J., "An abridged model of the
**     precession-nutation of the celestial pole", Cel.Mech.Dyn.Astron.
**     85, 37-49 (2003)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J., Astron.Astrophys. 282, 663-683 (1994)
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    _sofa.iauNut00b(date1,date2,byref(dpsi),byref(deps))
    return dpsi.value, deps.value


#iauNut06a
try:
    _sofa.iauNut06a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        ]
except  AttributeError:
    pass
def nut06a(date1,date2):
    """
    IAU 2000A nutation with adjustments to match the IAU 2006 precession.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsi : float
        nutation, luni-solar + planetary (Note 2).
    deps :float
        nutation, luni-solar + planetary (Note 2).
/*
**  - - - - - - - - - -
**   i a u N u t 0 6 a
**  - - - - - - - - - -
**
**  IAU 2000A nutation with adjustments to match the IAU 2006
**  precession.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps     double   nutation, luni-solar + planetary (Note 2)
**
**  Status:  canonical model.
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components in longitude and obliquity are in radians
**     and with respect to the mean equinox and ecliptic of date,
**     IAU 2006 precession model (Hilton et al. 2006, Capitaine et al.
**     2005).
**
**  3) The function first computes the IAU 2000A nutation, then applies
**     adjustments for (i) the consequences of the change in obliquity
**     from the IAU 1980 ecliptic to the IAU 2006 ecliptic and (ii) the
**     secular variation in the Earth's dynamical form factor J2.
**
**  4) The present function provides classical nutation, complementing
**     the IAU 2000 frame bias and IAU 2006 precession.  It delivers a
**     pole which is at current epochs accurate to a few tens of
**     microarcseconds, apart from the free core nutation.
**
**  Called:
**     iauNut00a    nutation, IAU 2000A
**
**  References:
**
**     Chapront, J., Chapront-Touze, M. & Francou, G. 2002,
**     Astron.Astrophys. 387, 700
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
**     107, B4.  The MHB_2000 code itself was obtained on 9th September
**     2002 from ftp//maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    _sofa.iauNut06a(date1,date2,byref(dpsi),byref(deps))
    return dpsi.value, deps.value


#iauNut80
try:
    _sofa.iauNut80.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        ]
except  AttributeError:
    pass
def nut80(date1,date2):
    """
    
Nutation, IAU 1980 model.
    Parameters
    ----------
    date1 : float
        DESCRIPTION.
    date2 : float
        DESCRIPTION.

    Returns
    -------
    dpsi : float
        nutation in longitude (radians).
    deps :float
        nutation in obliquity (radians).
/*
**  - - - - - - - - -
**   i a u N u t 8 0
**  - - - - - - - - -
**
**  Nutation, IAU 1980 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi          double    nutation in longitude (radians)
**     deps          double    nutation in obliquity (radians)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The nutation components are with respect to the ecliptic of
**     date.
**
**  Called:
**     iauAnpm      normalize angle into range +/- pi
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.222 (p111).
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    _sofa.iauNut80(date1,date2,byref(dpsi),byref(deps))
    return dpsi.value, deps.value


#iauNutm80
try:
    _sofa.iauNutm80.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape = (3,3), dtype = float) #rmatn[3][3]
        ]
except  AttributeError:
    pass
def nutm80(date1,date2):
    """
    Form the matrix of nutation for a given date, IAU 1980 model.

    Parameters
    ----------
    date1 : float
        TDB date (Note 1).
    date2 : float
        TDB date (Note 1).

    Returns
    -------
    rmatn : double[3][3]
        nutation matrix.
/*
**  - - - - - - - - - -
**   i a u N u t m 8 0
**  - - - - - - - - - -
**
**  Form the matrix of nutation for a given date, IAU 1980 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2    double          TDB date (Note 1)
**
**  Returned:
**     rmatn          double[3][3]    nutation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(true) = rmatn * V(mean),
**     where the p-vector V(true) is with respect to the true
**     equatorial triad of date and the p-vector V(mean) is with
**     respect to the mean equatorial triad of date.
**
**  Called:
**     iauNut80     nutation, IAU 1980
**     iauObl80     mean obliquity, IAU 1980
**     iauNumat     form nutation matrix
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatn = zeros(shape = (3,3), dtype = float)
    _sofa.iauNutm80(date1,date2,rmatn)
    return rmatn 


#iauObl06
try:
    _sofa.iauObl06.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauObl06.restype = c_double
except  AttributeError:
    pass
def obl06(date1,date2):
    """
    Mean obliquity of the ecliptic, IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    eps0 : float
        obliquity of the ecliptic (radians, Note 2).
/*
**  - - - - - - - - -
**   i a u O b l 0 6
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 2006 precession model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double   obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eps0 = _sofa.iauObl06(date1,date2)
    return eps0


#iauObl80
try:
    _sofa.iauObl80.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauObl80.restype = c_double
except  AttributeError:
    pass
def obl80(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    eps0 : float
        obliquity of the ecliptic (radians, Note 2).
/*
**  - - - - - - - - -
**   i a u O b l 8 0
**  - - - - - - - - -
**
**  Mean obliquity of the ecliptic, IAU 1980 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                   double    obliquity of the ecliptic (radians, Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The result is the angle between the ecliptic and mean equator of
**     date date1+date2.
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Expression 3.222-1 (p114).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eps0 = _sofa.iauObl80(date1,date2)
    return eps0

#iauPb06
try:
    _sofa.iauPb06.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#bzeta
        POINTER(c_double),#bz
        POINTER(c_double),#btheta
        ]
except  AttributeError:
    pass
def pb06(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    bzeta : float
        1st rotation: radians cw around z.
    bz : float
        3rd rotation: radians cw around z.
    btheta : float
        2nd rotation: radians ccw around y.
/*
**  - - - - - - - -
**   i a u P b 0 6
**  - - - - - - - -
**
**  This function forms three Euler angles which implement general
**  precession from epoch J2000.0, using the IAU 2006 model.  Frame
**  bias (the offset between ICRS and mean J2000.0) is included.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     bzeta        double   1st rotation: radians cw around z
**     bz           double   3rd rotation: radians cw around z
**     btheta       double   2nd rotation: radians ccw around y
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The traditional accumulated precession angles zeta_A, z_A,
**     theta_A cannot be obtained in the usual way, namely through
**     polynomial expressions, because of the frame bias.  The latter
**     means that two of the angles undergo rapid changes near this
**     date.  They are instead the results of decomposing the
**     precession-bias matrix obtained by using the Fukushima-Williams
**     method, which does not suffer from the problem.  The
**     decomposition returns values which can be used in the
**     conventional formulation and which include frame bias.
**
**  3) The three angles are returned in the conventional order, which
**     is not the same as the order of the corresponding Euler
**     rotations.  The precession-bias matrix is
**     R_3(-z) x R_2(+theta) x R_3(-zeta).
**
**  4) Should zeta_A, z_A, theta_A angles be required that do not
**     contain frame bias, they are available by calling the SOFA
**     function iauP06e.
**
**  Called:
**     iauPmat06    PB matrix, IAU 2006
**     iauRz        rotate around Z-axis
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    bzeta = c_double()
    bz = c_double()
    btheta = c_double()
    _sofa.iauPb06(date1,date2,byref(bzeta),byref(bz),byref(btheta))
    return bzeta.value, bz.value, btheta.value


#iauPfw06
try:
    _sofa.iauPfw06.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#gamb
        POINTER(c_double),#phib
        POINTER(c_double),#psib
        POINTER(c_double),#epsa
        ]
except  AttributeError:
    pass
def pfw06(date1,date2):
    """
    Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation)

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    float
        F-W angle gamma_bar (radians).
    float
        F-W angle phi_bar (radians).
    float
        F-W angle psi_bar (radians).
    float
        F-W angle epsilon_A (radians).
/*
**  - - - - - - - - -
**   i a u P f w 0 6
**  - - - - - - - - -
**
**  Precession angles, IAU 2006 (Fukushima-Williams 4-angle formulation).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     gamb         double   F-W angle gamma_bar (radians)
**     phib         double   F-W angle phi_bar (radians)
**     psib         double   F-W angle psi_bar (radians)
**     epsa         double   F-W angle epsilon_A (radians)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) Naming the following points:
**
**           e = J2000.0 ecliptic pole,
**           p = GCRS pole,
**           E = mean ecliptic pole of date,
**     and   P = mean pole of date,
**
**     the four Fukushima-Williams angles are as follows:
**
**        gamb = gamma_bar = epE
**        phib = phi_bar = pE
**        psib = psi_bar = pEP
**        epsa = epsilon_A = EP
**
**  3) The matrix representing the combined effects of frame bias and
**     precession is:
**
**        PxB = R_1(-epsa).R_3(-psib).R_1(phib).R_3(gamb)
**
**  4) The matrix representing the combined effects of frame bias,
**     precession and nutation is simply:
**
**        NxPxB = R_1(-epsa-dE).R_3(-psib-dP).R_1(phib).R_3(gamb)
**
**     where dP and dE are the nutation components with respect to the
**     ecliptic of date.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Called:
**     iauObl06     mean obliquity, IAU 2006
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    gamb = c_double()
    phib = c_double()
    psib = c_double()
    epsa = c_double()
    _sofa.iauPfw06(date1,date2,byref(gamb),byref(phib),byref(psib),byref(epsa))
    return gamb.value, phib.value, psib.value, epsa.value


#iauPmat00
try:
    _sofa.iauPmat00.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float)
        ]
except  AttributeError:
    pass
def pmat00(date1,date2):
    """
    Precession matrix (including frame bias) from GCRS to a specified date, IAU 2000 model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rbp : double[3][3]
        bias-precession matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u P m a t 0 0
**  - - - - - - - - - -
**
**  Precession matrix (including frame bias) from GCRS to a specified
**  date, IAU 2000 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbp          double[3][3]    bias-precession matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
**     the p-vector V(GCRS) is with respect to the Geocentric Celestial
**     Reference System (IAU, 2000) and the p-vector V(date) is with
**     respect to the mean equatorial triad of the given date.
**
**  Called:
**     iauBp00      frame bias and precession matrices, IAU 2000
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rbp = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPmat00(date1,date2,rbp)
    return rbp


#iauPmat06
try:
    _sofa.iauPmat06.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float)
        ]
except  AttributeError:
    pass
def pmat06(date1,date2):
    """
    Precession matrix (including frame bias) from GCRS to a specified  date, IAU 2006 model.

    Parameters
    ----------
    date1 : double
        TT as a 2-part Julian Date (Note 1).
    date2 : double
        bias-precession matrix (Note 2).

    Returns
    -------
    rbp : double[3][3]
        bias-precession matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u P m a t 0 6
**  - - - - - - - - - -
**
**  Precession matrix (including frame bias) from GCRS to a specified
**  date, IAU 2006 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbp          double[3][3]    bias-precession matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rbp * V(GCRS), where
**     the p-vector V(GCRS) is with respect to the Geocentric Celestial
**     Reference System (IAU, 2000) and the p-vector V(date) is with
**     respect to the mean equatorial triad of the given date.
**
**  Called:
**     iauPfw06     bias-precession F-W angles, IAU 2006
**     iauFw2m      F-W angles to r-matrix
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rbp = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPmat06(date1,date2,rbp)
    return rbp

#iauPmat76
try:
    _sofa.iauPmat76.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float), #rmatp[3][3]
        ]
except  AttributeError:
    pass
def pmat76(date1,date2):
    """
    Precession matrix from J2000.0 to a specified date, IAU 1976 model.

    Parameters
    ----------
    date1 : float
        ending date, TT (Note 1).
    date2 : float
        ending date, TT (Note 1).

    Returns
    -------
    rmatp : double[3][3]
        precession matrix, J2000.0 -> date1+date2.
/*
**  - - - - - - - - - -
**   i a u P m a t 7 6
**  - - - - - - - - - -
**
**  Precession matrix from J2000.0 to a specified date, IAU 1976 model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       ending date, TT (Note 1)
**
**  Returned:
**     rmatp       double[3][3] precession matrix, J2000.0 -> date1+date2
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = RMATP * V(J2000),
**     where the p-vector V(J2000) is with respect to the mean
**     equatorial triad of epoch J2000.0 and the p-vector V(date)
**     is with respect to the mean equatorial triad of the given
**     date.
**
**  3) Though the matrix method itself is rigorous, the precession
**     angles are expressed through canonical polynomials which are
**     valid only for a limited time span.  In addition, the IAU 1976
**     precession rate is known to be imperfect.  The absolute accuracy
**     of the present formulation is better than 0.1 arcsec from
**     1960AD to 2040AD, better than 1 arcsec from 1640AD to 2360AD,
**     and remains below 3 arcsec for the whole of the period
**     500BC to 3000AD.  The errors exceed 10 arcsec outside the
**     range 1200BC to 3900AD, exceed 100 arcsec outside 4200BC to
**     5600AD and exceed 1000 arcsec outside 6800BC to 8200AD.
**
**  Called:
**     iauPrec76    accumulated precession angles, IAU 1976
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRy        rotate around Y-axis
**     iauCr        copy r-matrix
**
**  References:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282.
**      equations (6) & (7), p283.
**
**     Kaplan,G.H., 1981. USNO circular no. 163, pA2.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatp = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPmat76(date1,date2,rmatp)
    return rmatp


#iauPn00
try:
    _sofa.iauPn00.argtypes = [
        c_double, #date1
        c_double, #date2
        c_double, #dpsi
        c_double, #deps
        POINTER(c_double),#epsa
        ndpointer(shape  = (3,3), dtype = float), #rb[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rn[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pn00(date1,date2,dpsi,dpes):
    """
    Precession-nutation, IAU 2000 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    dpsi : float
        nutation (Note 2).
    dpes : float
        nutation (Note 2).

    Returns
    -------
    epsa : float
        mean obliquity (Note 3).
    rb : double[3][3]
        frame bias matrix (Note 4).
    rp : double[3][3]
        precession matrix (Note 5).
    rbp : double[3][3]
        bias-precession matrix (Note 6).
    rn : double[3][3]
        nutation matrix (Note 7).
    rbpn : double[3][3]
        GCRS-to-true matrix (Note 8).
/*
**  - - - - - - - -
**   i a u P n 0 0
**  - - - - - - - -
**
**  Precession-nutation, IAU 2000 model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based
**  use indirectly.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**     dpsi,deps    double          nutation (Note 2)
**
**  Returned:
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Note 8)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The caller is responsible for providing the nutation components;
**     they are in longitude and obliquity, in radians and are with
**     respect to the equinox and ecliptic of date.  For high-accuracy
**     applications, free core nutation should be included as well as
**     any other relevant corrections to the position of the CIP.
**
**  3) The returned mean obliquity is consistent with the IAU 2000
**     precession-nutation models.
**
**  4) The matrix rb transforms vectors from GCRS to J2000.0 mean
**     equator and equinox by applying frame bias.
**
**  5) The matrix rp transforms vectors from J2000.0 mean equator and
**     equinox to mean equator and equinox of date by applying
**     precession.
**
**  6) The matrix rbp transforms vectors from GCRS to mean equator and
**     equinox of date by applying frame bias then precession.  It is
**     the product rp x rb.
**
**  7) The matrix rn transforms vectors from mean equator and equinox of
**     date to true equator and equinox of date by applying the nutation
**     (luni-solar + planetary).
**
**  8) The matrix rbpn transforms vectors from GCRS to true equator and
**     equinox of date.  It is the product rn x rbp, applying frame
**     bias, precession and nutation in that order.
**
**  9) It is permissible to re-use the same array in the returned
**     arguments.  The arrays are filled in the order given.
**
**  Called:
**     iauPr00      IAU 2000 precession adjustments
**     iauObl80     mean obliquity, IAU 1980
**     iauBp00      frame bias and precession matrices, IAU 2000
**     iauCr        copy r-matrix
**     iauNumat     form nutation matrix
**     iauRxr       product of two r-matrices
**
**  Reference:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    epsa = c_double()
    rb = zeros(shape = (3,3), dtype =  float)
    rp = zeros(shape = (3,3), dtype =  float)
    rbp = zeros(shape = (3,3), dtype =  float)
    rn = zeros(shape = (3,3), dtype =  float)
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPn00(date1,date2,dpsi,dpes,byref(epsa),rb,rp,rbp,rn,rbpn)
    return epsa.value, rb, rp, rbp, rn, rbpn


#iauPn00a
try:
    _sofa.iauPn00a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        POINTER(c_double),#epsa
        ndpointer(shape  = (3,3), dtype = float), #rb[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rn[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pn00a(date1,date2):
    """
    Precession-nutation, IAU 2000A model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsi : float
        nutation (Note 2).
    deps : float
        nutation (Note 2).
    epsa :float
        mean obliquity (Note 3).
    rb : double[3][3]
        frame bias matrix (Note 4).
    rp : double[3][3]
        precession matrix (Note 5).
    rbp : double[3][3]
        bias-precession matrix (Note 6).
    rn : double[3][3]
        nutation matrix (Note 7).
    rbpn : double[3][3]
        GCRS-to-true matrix (Notes 8,9).

    """
    dpsi = c_double()
    deps = c_double()
    epsa = c_double()
    rb = zeros(shape = (3,3), dtype =  float)
    rp = zeros(shape = (3,3), dtype =  float)
    rbp = zeros(shape = (3,3), dtype =  float)
    rn = zeros(shape = (3,3), dtype =  float)
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPn00a(date1,date2,byref(dpsi),byref(deps),byref(epsa),rb,rp,rbp,rn,rbpn)
    return dpsi.value,deps.value,epsa.value, rb, rp, rbp, rn, rbpn


#iauPn00b
try:
    _sofa.iauPn00b.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        POINTER(c_double),#epsa
        ndpointer(shape  = (3,3), dtype = float), #rb[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rn[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pn00b(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
         TT as a 2-part Julian Date (Note 1).
    date2 : float
         TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsi : float
        nutation (Note 2).
    deps : float
        nutation (Note 2).
    epsa : float
        mean obliquity (Note 3).
    rb : double[3][3]
        frame bias matrix (Note 4).
    rp : double[3][3]
        precession matrix (Note 5).
    rbp : double[3][3]
        bias-precession matrix (Note 6).
    rn : double[3][3]
        nutation matrix (Note 7).
    rbpn : double[3][3]
        GCRS-to-true matrix (Notes 8,9).
/*
**  - - - - - - - - -
**   i a u P n 0 0 b
**  - - - - - - - - -
**
**  Precession-nutation, IAU 2000B model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based
**  use indirectly.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps    double          nutation (Note 2)
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The nutation components (luni-solar + planetary, IAU 2000B) in
**      longitude and obliquity are in radians and with respect to the
**      equinox and ecliptic of date.  For more accurate results, but
**      at the cost of increased computation, use the iauPn00a function.
**      For the utmost accuracy, use the iauPn00  function, where the
**      nutation components are caller-specified.
**
**  3)  The mean obliquity is consistent with the IAU 2000 precession.
**
**  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
**      equator and equinox by applying frame bias.
**
**  5)  The matrix rp transforms vectors from J2000.0 mean equator and
**      equinox to mean equator and equinox of date by applying
**      precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean equator and
**      equinox of date by applying frame bias then precession.  It is
**      the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean equator and equinox
**      of date to true equator and equinox of date by applying the
**      nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true equator and
**      equinox of date.  It is the product rn x rbp, applying frame
**      bias, precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2000B Celestial Intermediate
**      Pole are elements (3,1-3) of the GCRS-to-true matrix,
**      i.e. rbpn[2][0-2].
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     iauNut00b    nutation, IAU 2000B
**     iauPn00      bias/precession/nutation results, IAU 2000
**
**  Reference:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003).
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**  This revision:  2013 November 13
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    epsa = c_double()
    rb = zeros(shape = (3,3), dtype =  float)
    rp = zeros(shape = (3,3), dtype =  float)
    rbp = zeros(shape = (3,3), dtype =  float)
    rn = zeros(shape = (3,3), dtype =  float)
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPn00b(date1,date2,byref(dpsi),byref(deps),byref(epsa),rb,rp,rbp,rn,rbpn)
    return dpsi.value,deps.value,epsa.value, rb, rp, rbp, rn, rbpn


#iauPn06
try:
    _sofa.iauPn06.argtypes = [
        c_double, #date1
        c_double, #date2
        c_double, #dpsi
        c_double, #deps
        POINTER(c_double),#epsa
        ndpointer(shape  = (3,3), dtype = float), #rb[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rn[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pn06(date1,date2,dpsi,dpes):
    """
    Precession-nutation, IAU 2006 model:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    dpsi : float
        nutation (Note 2).
    dpes : float
        nutation (Note 2).

    Returns
    -------
    epsa : float
        mean obliquity (Note 3).
    rb : double[3][3]
        frame bias matrix (Note 4).
    rp : double[3][3]
        precession matrix (Note 5).
    rbp : double[3][3]
        DESCRbias-precession matrix (Note 6)IPTION.
    rn : double[3][3]
        nutation matrix (Note 7).
    rbpn : double[3][3]
        GCRS-to-true matrix (Note 8).
/*
**  - - - - - - - -
**   i a u P n 0 6
**  - - - - - - - -
**
**  Precession-nutation, IAU 2006 model:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based use
**  indirectly.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**     dpsi,deps    double          nutation (Note 2)
**
**  Returned:
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Note 8)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The caller is responsible for providing the nutation components;
**      they are in longitude and obliquity, in radians and are with
**      respect to the equinox and ecliptic of date.  For high-accuracy
**      applications, free core nutation should be included as well as
**      any other relevant corrections to the position of the CIP.
**
**  3)  The returned mean obliquity is consistent with the IAU 2006
**      precession.
**
**  4)  The matrix rb transforms vectors from GCRS to J2000.0 mean
**      equator and equinox by applying frame bias.
**
**  5)  The matrix rp transforms vectors from J2000.0 mean equator and
**      equinox to mean equator and equinox of date by applying
**      precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean equator and
**      equinox of date by applying frame bias then precession.  It is
**      the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean equator and equinox
**      of date to true equator and equinox of date by applying the
**      nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true equator and
**      equinox of date.  It is the product rn x rbp, applying frame
**      bias, precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the Celestial Intermediate Pole are
**      elements (3,1-3) of the GCRS-to-true matrix, i.e. rbpn[2][0-2].
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     iauPfw06     bias-precession F-W angles, IAU 2006
**     iauFw2m      F-W angles to r-matrix
**     iauCr        copy r-matrix
**     iauTr        transpose r-matrix
**     iauRxr       product of two r-matrices
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 November 14
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    epsa = c_double()
    rb = zeros(shape = (3,3), dtype =  float)
    rp = zeros(shape = (3,3), dtype =  float)
    rbp = zeros(shape = (3,3), dtype =  float)
    rn = zeros(shape = (3,3), dtype =  float)
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPn06(date1,date2,dpsi,dpes,byref(epsa),rb,rp,rbp,rn,rbpn)
    return epsa.value, rb, rp, rbp, rn, rbpn


#iauPn06a
try:
    _sofa.iauPn06a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsi
        POINTER(c_double), #deps
        POINTER(c_double),#epsa
        ndpointer(shape  = (3,3), dtype = float), #rb[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbp[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rn[3][3]
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pn06a(date1,date2):
    """
    Precession-nutation, IAU 2006/2000A models:  a multi-purpose function, supporting classical (equinox-based) use directly and CIO-based use indirectly.

    Parameters
    ----------
    date1 : float
        DESCRIPTION.
    date2 : float
        DESCRIPTION.

    Returns
    -------
    dpsi : float
        nutation (Note 2).
    deps : float
        nutation (Note 2).
    epsa : float
        mean obliquity (Note 3).
    rb : TYPE
        frame bias matrix (Note 4).
    rp : TYPE
        precession matrix (Note 5).
    rbp : TYPE
        bias-precession matrix (Note 6).
    rn : TYPE
        nutation matrix (Note 7).
    rbpn : TYPE
        GCRS-to-true matrix (Notes 8,9).
/*
**  - - - - - - - - -
**   i a u P n 0 6 a
**  - - - - - - - - -
**
**  Precession-nutation, IAU 2006/2000A models:  a multi-purpose function,
**  supporting classical (equinox-based) use directly and CIO-based use
**  indirectly.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double          TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsi,deps    double          nutation (Note 2)
**     epsa         double          mean obliquity (Note 3)
**     rb           double[3][3]    frame bias matrix (Note 4)
**     rp           double[3][3]    precession matrix (Note 5)
**     rbp          double[3][3]    bias-precession matrix (Note 6)
**     rn           double[3][3]    nutation matrix (Note 7)
**     rbpn         double[3][3]    GCRS-to-true matrix (Notes 8,9)
**
**  Notes:
**
**  1)  The TT date date1+date2 is a Julian Date, apportioned in any
**      convenient way between the two arguments.  For example,
**      JD(TT)=2450123.7 could be expressed in any of these ways,
**      among others:
**
**             date1          date2
**
**          2450123.7           0.0       (JD method)
**          2451545.0       -1421.3       (J2000 method)
**          2400000.5       50123.2       (MJD method)
**          2450123.5           0.2       (date & time method)
**
**      The JD method is the most natural and convenient to use in
**      cases where the loss of several decimal digits of resolution
**      is acceptable.  The J2000 method is best matched to the way
**      the argument is handled internally and will deliver the
**      optimum resolution.  The MJD method and the date & time methods
**      are both good compromises between resolution and convenience.
**
**  2)  The nutation components (luni-solar + planetary, IAU 2000A) in
**      longitude and obliquity are in radians and with respect to the
**      equinox and ecliptic of date.  Free core nutation is omitted;
**      for the utmost accuracy, use the iauPn06 function, where the
**      nutation components are caller-specified.
**
**  3)  The mean obliquity is consistent with the IAU 2006 precession.
**
**  4)  The matrix rb transforms vectors from GCRS to mean J2000.0 by
**      applying frame bias.
**
**  5)  The matrix rp transforms vectors from mean J2000.0 to mean of
**      date by applying precession.
**
**  6)  The matrix rbp transforms vectors from GCRS to mean of date by
**      applying frame bias then precession.  It is the product rp x rb.
**
**  7)  The matrix rn transforms vectors from mean of date to true of
**      date by applying the nutation (luni-solar + planetary).
**
**  8)  The matrix rbpn transforms vectors from GCRS to true of date
**      (CIP/equinox).  It is the product rn x rbp, applying frame bias,
**      precession and nutation in that order.
**
**  9)  The X,Y,Z coordinates of the IAU 2006/2000A Celestial
**      Intermediate Pole are elements (3,1-3) of the GCRS-to-true
**      matrix, i.e. rbpn[2][0-2].
**
**  10) It is permissible to re-use the same array in the returned
**      arguments.  The arrays are filled in the stated order.
**
**  Called:
**     iauNut06a    nutation, IAU 2006/2000A
**     iauPn06      bias/precession/nutation results, IAU 2006
**
**  Reference:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**  This revision:  2013 November 13
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsi = c_double()
    deps = c_double()
    epsa = c_double()
    rb = zeros(shape = (3,3), dtype =  float)
    rp = zeros(shape = (3,3), dtype =  float)
    rbp = zeros(shape = (3,3), dtype =  float)
    rn = zeros(shape = (3,3), dtype =  float)
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPn06a(date1,date2,byref(dpsi),byref(deps),byref(epsa),rb,rp,rbp,rn,rbpn)
    return dpsi.value,deps.value,epsa.value, rb, rp, rbp, rn, rbpn


#iauPnm00a
try:
    _sofa.iauPnm00a.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pnm00a(date1,date2):
    """
    Form the matrix of precession-nutation for a given date (including  frame bias), equinox-based, IAU 2000A model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rbpn : double[3][3]
        classical NPB matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u P n m 0 0 a
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), equinox-based, IAU 2000A model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double     TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbpn         double[3][3]    classical NPB matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  3) A faster, but slightly less accurate result (about 1 mas), can be
**     obtained by using instead the iauPnm00b function.
**
**  Called:
**     iauPn00a     bias/precession/nutation, IAU 2000A
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPnm00a(date1,date2,rbpn)
    return rbpn


#iauPnm00b
try:
    _sofa.iauPnm00b.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float), #rbpn[3][3]
        ]
except  AttributeError:
    pass
def pnm00b(date1,date2):
    """
Form the matrix of precession-nutation for a given date (including frame bias), equinox-based, IAU 2000B model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rbpn : double[3][3]
        bias-precession-nutation matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u P n m 0 0 b
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), equinox-based, IAU 2000B model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rbpn        double[3][3] bias-precession-nutation matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rbpn * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  3) The present function is faster, but slightly less accurate (about
**     1 mas), than the iauPnm00a function.
**
**  Called:
**     iauPn00b     bias/precession/nutation, IAU 2000B
**
**  Reference:
**
**     IAU: Trans. International Astronomical Union, Vol. XXIVB;  Proc.
**     24th General Assembly, Manchester, UK.  Resolutions B1.3, B1.6.
**     (2000)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rbpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPnm00b(date1,date2,rbpn)
    return rbpn


#iauPnm06a
try:
    _sofa.iauPnm06a.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float), #rnpb[3][3]
        ]
except  AttributeError:
    pass
def pnm06a(date1,date2):
    """
    Form the matrix of precession-nutation for a given date (including frame bias), IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    rnpb : double[3][3]
        bias-precession-nutation matrix (Note 2).
/*
**  - - - - - - - - - -
**   i a u P n m 0 6 a
**  - - - - - - - - - -
**
**  Form the matrix of precession-nutation for a given date (including
**  frame bias), IAU 2006 precession and IAU 2000A nutation models.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double       TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     rnpb        double[3][3] bias-precession-nutation matrix (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rnpb * V(GCRS), where
**     the p-vector V(date) is with respect to the true equatorial triad
**     of date date1+date2 and the p-vector V(GCRS) is with respect to
**     the Geocentric Celestial Reference System (IAU, 2000).
**
**  Called:
**     iauPfw06     bias-precession F-W angles, IAU 2006
**     iauNut06a    nutation, IAU 2006/2000A
**     iauFw2m      F-W angles to r-matrix
**
**  Reference:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rnpb = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPnm06a(date1,date2,rnpb)
    return rnpb


#iauPnm80
try:
    _sofa.iauPnm80.argtypes = [
        c_double, #date1
        c_double, #date2
        ndpointer(shape  = (3,3), dtype = float), #rmatpn[3][3]
        ]
except  AttributeError:
    pass
def pnm80(date1,date2):
    """
    Form the matrix of precession/nutation for a given date, IAU 1976  precession model, IAU 1980 nutation model.

    Parameters
    ----------
    date1 : float
        TDB date (Note 1).
    date2 : float
        TDB date (Note 1).

    Returns
    -------
    rmatpn : double[3][3]
        combined precession/nutation matrix.
/*
**  - - - - - - - - -
**   i a u P n m 8 0
**  - - - - - - - - -
**
**  Form the matrix of precession/nutation for a given date, IAU 1976
**  precession model, IAU 1980 nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2    double         TDB date (Note 1)
**
**  Returned:
**     rmatpn         double[3][3]   combined precession/nutation matrix
**
**  Notes:
**
**  1) The TDB date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TDB)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The matrix operates in the sense V(date) = rmatpn * V(J2000),
**     where the p-vector V(date) is with respect to the true equatorial
**     triad of date date1+date2 and the p-vector V(J2000) is with
**     respect to the mean equatorial triad of epoch J2000.0.
**
**  Called:
**     iauPmat76    precession matrix, IAU 1976
**     iauNutm80    nutation matrix, IAU 1980
**     iauRxr       product of two r-matrices
**
**  Reference:
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 3.3 (p145).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rmatpn = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPnm80(date1,date2,rmatpn)
    return rmatpn


#iauP06e
try:
    _sofa.iauP06e.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#eps0
        POINTER(c_double),#psia
        POINTER(c_double),#oma
        POINTER(c_double),#bpa
        POINTER(c_double),#bqa
        POINTER(c_double),#pia
        POINTER(c_double),#bpia
        POINTER(c_double),#epsa
        POINTER(c_double),#chia
        POINTER(c_double),#za
        POINTER(c_double),#zetaa
        POINTER(c_double),#thetaa
        POINTER(c_double),#pa
        POINTER(c_double),#gam
        POINTER(c_double),#phi
        POINTER(c_double),#psi
        ]
except  AttributeError:
    pass
def p06e(date1,date2):
    """
    Precession angles, IAU 2006, equinox based.

    Parameters
    ----------
    date1 : flaot
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    eps0 : float
        epsilon_0.
    psia : float
        psi_A.
    oma : float
        omega_A.
    bpa : float
        P_A.
    bqa : float
        Q_A.
    pia : float
        pi_A.
    bpia : float
        Pi_A.
    epsa : float
        obliquity epsilon_A.
    chia : float
        chi_A.
    za : float
        z_A.
    zetaa : float
        zeta_A.
    thetaa : float
        theta_A.
    pa : float
        p_A.
    gam : float
        F-W angle gamma_J2000.
    phi : float
        F-W angle phi_J2000.
    psi : float
        F-W angle psi_J2000.
/*
**  - - - - - - - -
**   i a u P 0 6 e
**  - - - - - - - -
**
**  Precession angles, IAU 2006, equinox based.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical models.
**
**  Given:
**     date1,date2   double   TT as a 2-part Julian Date (Note 1)
**
**  Returned (see Note 2):
**     eps0          double   epsilon_0
**     psia          double   psi_A
**     oma           double   omega_A
**     bpa           double   P_A
**     bqa           double   Q_A
**     pia           double   pi_A
**     bpia          double   Pi_A
**     epsa          double   obliquity epsilon_A
**     chia          double   chi_A
**     za            double   z_A
**     zetaa         double   zeta_A
**     thetaa        double   theta_A
**     pa            double   p_A
**     gam           double   F-W angle gamma_J2000
**     phi           double   F-W angle phi_J2000
**     psi           double   F-W angle psi_J2000
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) This function returns the set of equinox based angles for the
**     Capitaine et al. "P03" precession theory, adopted by the IAU in
**     2006.  The angles are set out in Table 1 of Hilton et al. (2006):
**
**     eps0   epsilon_0   obliquity at J2000.0
**     psia   psi_A       luni-solar precession
**     oma    omega_A     inclination of equator wrt J2000.0 ecliptic
**     bpa    P_A         ecliptic pole x, J2000.0 ecliptic triad
**     bqa    Q_A         ecliptic pole -y, J2000.0 ecliptic triad
**     pia    pi_A        angle between moving and J2000.0 ecliptics
**     bpia   Pi_A        longitude of ascending node of the ecliptic
**     epsa   epsilon_A   obliquity of the ecliptic
**     chia   chi_A       planetary precession
**     za     z_A         equatorial precession: -3rd 323 Euler angle
**     zetaa  zeta_A      equatorial precession: -1st 323 Euler angle
**     thetaa theta_A     equatorial precession: 2nd 323 Euler angle
**     pa     p_A         general precession
**     gam    gamma_J2000 J2000.0 RA difference of ecliptic poles
**     phi    phi_J2000   J2000.0 codeclination of ecliptic pole
**     psi    psi_J2000   longitude difference of equator poles, J2000.0
**
**     The returned values are all radians.
**
**  3) Hilton et al. (2006) Table 1 also contains angles that depend on
**     models distinct from the P03 precession theory itself, namely the
**     IAU 2000A frame bias and nutation.  The quoted polynomials are
**     used in other SOFA functions:
**
**     . iauXy06  contains the polynomial parts of the X and Y series.
**
**     . iauS06  contains the polynomial part of the s+XY/2 series.
**
**     . iauPfw06  implements the series for the Fukushima-Williams
**       angles that are with respect to the GCRS pole (i.e. the variants
**       that include frame bias).
**
**  4) The IAU resolution stipulated that the choice of parameterization
**     was left to the user, and so an IAU compliant precession
**     implementation can be constructed using various combinations of
**     the angles returned by the present function.
**
**  5) The parameterization used by SOFA is the version of the Fukushima-
**     Williams angles that refers directly to the GCRS pole.  These
**     angles may be calculated by calling the function iauPfw06.  SOFA
**     also supports the direct computation of the CIP GCRS X,Y by
**     series, available by calling iauXy06.
**
**  6) The agreement between the different parameterizations is at the
**     1 microarcsecond level in the present era.
**
**  7) When constructing a precession formulation that refers to the GCRS
**     pole rather than the dynamical pole, it may (depending on the
**     choice of angles) be necessary to introduce the frame bias
**     explicitly.
**
**  8) It is permissible to re-use the same variable in the returned
**     arguments.  The quantities are stored in the stated order.
**
**  Reference:
**
**     Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351
**
**  Called:
**     iauObl06     mean obliquity, IAU 2006
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    eps0 = c_double()
    psia = c_double()
    oma = c_double()
    bpa = c_double()
    bqa = c_double()
    pia = c_double()
    bpia = c_double()
    epsa = c_double()
    chia = c_double()
    za = c_double()
    zetaa = c_double()
    thetaa = c_double()
    pa = c_double()
    gam = c_double()
    phi = c_double()
    psi = c_double()
    _sofa.iauP06e(date1,date2,byref(eps0),byref(psia),byref(oma),byref(bpa),byref(bqa),byref(pia),byref(bpia),byref(epsa),byref(chia),byref(za),byref(zetaa),byref(thetaa),byref(pa),byref(gam),byref(phi),byref(psi))
    return eps0.value,psia.value,oma.value,bpa.value,bqa.value,pia.value,bpia.value,epsa.value,chia.value,za.value,zetaa.value,thetaa.value,pa.value,gam.value,phi.value,psi.value


#iauPom00
try:
    _sofa.iauPom00.argtypes = [
        c_double, #xp
        c_double, #yp
        c_double, #sp
        ndpointer(shape  = (3,3), dtype = float), #rpom[3][3]
        ]
except  AttributeError:
    pass
def pom00(xp,yp,sp):
    """
    

    Parameters
    ----------
    xp : float
        coordinates of the pole (radians, Note 1).
    yp : float
        coordinates of the pole (radians, Note 1).
    sp : flaot
        the TIO locator s' (radians, Note 2).

    Returns
    -------
    rpom : double[3][3]
        polar-motion matrix (Note 3).
/*
**  - - - - - - - - - -
**   i a u P o m 0 0
**  - - - - - - - - - -
**
**  Form the matrix of polar motion for a given date, IAU 2000.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xp,yp    double    coordinates of the pole (radians, Note 1)
**     sp       double    the TIO locator s' (radians, Note 2)
**
**  Returned:
**     rpom     double[3][3]   polar-motion matrix (Note 3)
**
**  Notes:
**
**  1) The arguments xp and yp are the coordinates (in radians) of the
**     Celestial Intermediate Pole with respect to the International
**     Terrestrial Reference System (see IERS Conventions 2003),
**     measured along the meridians to 0 and 90 deg west respectively.
**
**  2) The argument sp is the TIO locator s', in radians, which
**     positions the Terrestrial Intermediate Origin on the equator.  It
**     is obtained from polar motion observations by numerical
**     integration, and so is in essence unpredictable.  However, it is
**     dominated by a secular drift of about 47 microarcseconds per
**     century, and so can be taken into account by using s' = -47*t,
**     where t is centuries since J2000.0.  The function iauSp00
**     implements this approximation.
**
**  3) The matrix operates in the sense V(TRS) = rpom * V(CIP), meaning
**     that it is the final rotation when computing the pointing
**     direction to a celestial source.
**
**  Called:
**     iauIr        initialize r-matrix to identity
**     iauRz        rotate around Z-axis
**     iauRy        rotate around Y-axis
**     iauRx        rotate around X-axis
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rpom = zeros(shape = (3,3), dtype =  float)
    _sofa.iauPom00(xp,yp,sp,rpom)
    return rpom


#iauPr00
try:
    _sofa.iauPr00.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double), #dpsipr
        POINTER(c_double), #depspr
        ]
except  AttributeError:
    pass
def pr00(date1,date2):
    """
    Precession-rate part of the IAU 2000 precession-nutation models  (part of MHB2000).

    Parameters
    ----------
    date1 : float
        double  TT as a 2-part Julian Date (Note 1).
    date2 : flaot
        double  TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    dpsipr : float
        precession corrections (Notes 2,3).
    depspr : float
        precession corrections (Notes 2,3).
/*
**  - - - - - - - -
**   i a u P r 0 0
**  - - - - - - - -
**
**  Precession-rate part of the IAU 2000 precession-nutation models
**  (part of MHB2000).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2    double  TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     dpsipr,depspr  double  precession corrections (Notes 2,3)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The precession adjustments are expressed as "nutation
**     components", corrections in longitude and obliquity with respect
**     to the J2000.0 equinox and ecliptic.
**
**  3) Although the precession adjustments are stated to be with respect
**     to Lieske et al. (1977), the MHB2000 model does not specify which
**     set of Euler angles are to be used and how the adjustments are to
**     be applied.  The most literal and straightforward procedure is to
**     adopt the 4-rotation epsilon_0, psi_A, omega_A, xi_A option, and
**     to add dpsipr to psi_A and depspr to both omega_A and eps_A.
**
**  4) This is an implementation of one aspect of the IAU 2000A nutation
**     model, formally adopted by the IAU General Assembly in 2000,
**     namely MHB2000 (Mathews et al. 2002).
**
**  References:
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B., "Expressions
**     for the precession quantities based upon the IAU (1976) System of
**     Astronomical Constants", Astron.Astrophys., 58, 1-16 (1977)
**
**     Mathews, P.M., Herring, T.A., Buffet, B.A., "Modeling of nutation
**     and precession   New nutation series for nonrigid Earth and
**     insights into the Earth's interior", J.Geophys.Res., 107, B4,
**     2002.  The MHB2000 code itself was obtained on 9th September 2002
**     from ftp://maia.usno.navy.mil/conv2000/chapter5/IAU2000A.
**
**     Wallace, P.T., "Software for Implementing the IAU 2000
**     Resolutions", in IERS Workshop 5.1 (2002).
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dpsipr = c_double()
    depspr = c_double()
    _sofa.iauPr00(date1,date2,byref(dpsipr),byref(depspr))
    return dpsipr.value, depspr.value


#iauPrec76
try:
    _sofa.iauPrec76.argtypes = [
        c_double, #date01
        c_double, #date02
        c_double, #date11
        c_double, #date12
        POINTER(c_double), #zeta
        POINTER(c_double), #z
        POINTER(c_double), #theta
        ]
except  AttributeError:
    pass
def prec76(date01,date02,date11,date12):
    """
    This function forms the three Euler angles which implement general precession between two dates, using the IAU 1976 model (as for the FK5 catalog).

    Parameters
    ----------
    date01 : float
        TDB starting date (Note 1).
    date02 : float
        DESCRITDB starting date (Note 1)PTION.
    date11 : flaot
        TDB ending date (Note 1).
    date12 : flaot
        TDB ending date (Note 1).

    Returns
    -------
    float
        1st rotation: radians cw around z.
    float
        3rd rotation: radians cw around z.
    float
        2nd rotation: radians ccw around y.
/*
**  - - - - - - - - - -
**   i a u P r e c 7 6
**  - - - - - - - - - -
**
**  IAU 1976 precession model.
**
**  This function forms the three Euler angles which implement general
**  precession between two dates, using the IAU 1976 model (as for the
**  FK5 catalog).
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date01,date02   double    TDB starting date (Note 1)
**     date11,date12   double    TDB ending date (Note 1)
**
**  Returned:
**     zeta            double    1st rotation: radians cw around z
**     z               double    3rd rotation: radians cw around z
**     theta           double    2nd rotation: radians ccw around y
**
**  Notes:
**
**  1) The dates date01+date02 and date11+date12 are Julian Dates,
**     apportioned in any convenient way between the arguments daten1
**     and daten2.  For example, JD(TDB)=2450123.7 could be expressed in
**     any of these ways, among others:
**
**           daten1        daten2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in cases
**     where the loss of several decimal digits of resolution is
**     acceptable.  The J2000 method is best matched to the way the
**     argument is handled internally and will deliver the optimum
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**     The two dates may be expressed using different methods, but at
**     the risk of losing some resolution.
**
**  2) The accumulated precession angles zeta, z, theta are expressed
**     through canonical polynomials which are valid only for a limited
**     time span.  In addition, the IAU 1976 precession rate is known to
**     be imperfect.  The absolute accuracy of the present formulation
**     is better than 0.1 arcsec from 1960AD to 2040AD, better than
**     1 arcsec from 1640AD to 2360AD, and remains below 3 arcsec for
**     the whole of the period 500BC to 3000AD.  The errors exceed
**     10 arcsec outside the range 1200BC to 3900AD, exceed 100 arcsec
**     outside 4200BC to 5600AD and exceed 1000 arcsec outside 6800BC to
**     8200AD.
**
**  3) The three angles are returned in the conventional order, which
**     is not the same as the order of the corresponding Euler
**     rotations.  The precession matrix is
**     R_3(-z) x R_2(+theta) x R_3(-zeta).
**
**  Reference:
**
**     Lieske, J.H., 1979, Astron.Astrophys. 73, 282, equations
**     (6) & (7), p283.
**
**  This revision:  2013 November 19
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    zeta = c_double()
    z = c_double()
    theta = c_double()
    _sofa.iauPrec76(date01,date02,date11,date12,byref(zeta),byref(z),byref(theta))
    return zeta.value, z.value, theta.value


#iauS00
try:
    _sofa.iauS00.argtypes = [
        c_double, #date1
        c_double, #date2
        c_double, #x
        c_double, #y
        ]
    _sofa.iauS00.restype = c_double
except  AttributeError:
    pass
def s00(date1,date2,x,y):
    """
    The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y  coordinates.  Compatible with IAU 2000A precession-nutation.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    x : float
        CIP coordinates (Note 3).
    y : float
        CIP coordinates (Note 3).

    Returns
    -------
    s : float
        the CIO locator s in radians (Note 2).
/*
**  - - - - - - -
**   i a u S 0 0
**  - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
**  coordinates.  Compatible with IAU 2000A precession-nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**     x,y           double    CIP coordinates (Note 3)
**
**  Returned (function value):
**                   double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems:  the two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The quantity s remains below 0.1 arcsecond
**     throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  This
**     function requires X,Y to be supplied by the caller, who is
**     responsible for providing values that are consistent with the
**     supplied date.
**
**  4) The model is consistent with the IAU 2000A precession-nutation.
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFalp03    mean anomaly of the Sun
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFad03     mean elongation of the Moon from the Sun
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauS00(date1,date2,x,y)
    return s


#iauS00a
try:
    _sofa.iauS00a.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauS00a.restype = c_double
except  AttributeError:
    pass
def s00a(date1,date2):
    """
    The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2000A precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    s : float
        the CIO locator s in radians (Note 2).
/*
**  - - - - - - - -
**   i a u S 0 0 a
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2000A
**  precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  The present
**     function uses the full IAU 2000A nutation model when predicting
**     the CIP position.  Faster results, with no significant loss of
**     accuracy, can be obtained via the function iauS00b, which uses
**     instead the IAU 2000B truncated model.
**
**  Called:
**     iauPnm00a    classical NPB matrix, IAU 2000A
**     iauBnp2xy    extract CIP X,Y from the BPN matrix
**     iauS00       the CIO locator s, given X,Y, IAU 2000A
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauS00a(date1,date2)
    return s


#iauS00b
try:
    _sofa.iauS00b.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauS00b.restype = c_double
except  AttributeError:
    pass
def s00b(date1,date2):
    """
    The CIO locator s, positioning the Celestial Intermediate Origin on  the equator of the Celestial Intermediate Pole, using the IAU 2000B precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    s : float
        the CIO locator s in radians (Note 2).
/*
**  - - - - - - - -
**   i a u S 0 0 b
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2000B
**  precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  The present
**     function uses the IAU 2000B truncated nutation model when
**     predicting the CIP position.  The function iauS00a uses instead
**     the full IAU 2000A model, but with no significant increase in
**     accuracy and at some cost in speed.
**
**  Called:
**     iauPnm00b    classical NPB matrix, IAU 2000B
**     iauBnp2xy    extract CIP X,Y from the BPN matrix
**     iauS00       the CIO locator s, given X,Y, IAU 2000A
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauS00b(date1,date2)
    return s


#iauS06
try:
    _sofa.iauS06.argtypes = [
        c_double, #date1
        c_double, #date2
        c_double, #x
        c_double, #y
        ]
    _sofa.iauS06.restype = c_double
except  AttributeError:
    pass
def s06(date1,date2,x,y):
    """
   The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, given the CIP's X,Y coordinates.  Compatible with IAU 2006/2000A precession-nutation.


    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).
    x : float
        CIP coordinates (Note 3).
    y : float
        CIP coordinates (Note 3).

    Returns
    -------
    s : float
        the CIO locator s in radians (Note 2).
/*
**  - - - - - - -
**   i a u S 0 6
**  - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, given the CIP's X,Y
**  coordinates.  Compatible with IAU 2006/2000A precession-nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2   double    TT as a 2-part Julian Date (Note 1)
**     x,y           double    CIP coordinates (Note 3)
**
**  Returned (function value):
**                   double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems:  the two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The quantity s remains below 0.1 arcsecond
**     throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series
**     is more compact than a direct series for s would be.  This
**     function requires X,Y to be supplied by the caller, who is
**     responsible for providing values that are consistent with the
**     supplied date.
**
**  4) The model is consistent with the "P03" precession (Capitaine et
**     al. 2003), adopted by IAU 2006 Resolution 1, 2006, and the
**     IAU 2000A nutation (with P03 adjustments).
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFalp03    mean anomaly of the Sun
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFad03     mean elongation of the Moon from the Sun
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2003, Astron.
**     Astrophys. 432, 355
**
**     McCarthy, D.D., Petit, G. (eds.) 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauS06(date1,date2,x,y)
    return s


#iauS06a
try:
    _sofa.iauS06a.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauS06a.restype = c_double
except  AttributeError:
    pass
def s06a(date1,date2):
    """
The CIO locator s, positioning the Celestial Intermediate Origin on the equator of the Celestial Intermediate Pole, using the IAU 2006 precession and IAU 2000A nutation models.
    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    s : float
        the CIO locator s in radians (Note 2).
/*
**  - - - - - - - -
**   i a u S 0 6 a
**  - - - - - - - -
**
**  The CIO locator s, positioning the Celestial Intermediate Origin on
**  the equator of the Celestial Intermediate Pole, using the IAU 2006
**  precession and IAU 2000A nutation models.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the CIO locator s in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The CIO locator s is the difference between the right ascensions
**     of the same point in two systems.  The two systems are the GCRS
**     and the CIP,CIO, and the point is the ascending node of the
**     CIP equator.  The CIO locator s remains a small fraction of
**     1 arcsecond throughout 1900-2100.
**
**  3) The series used to compute s is in fact for s+XY/2, where X and Y
**     are the x and y components of the CIP unit vector;  this series is
**     more compact than a direct series for s would be.  The present
**     function uses the full IAU 2000A nutation model when predicting
**     the CIP position.
**
**  Called:
**     iauPnm06a    classical NPB matrix, IAU 2006/2000A
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS06       the CIO locator s, given X,Y, IAU 2006
**
**  References:
**
**     Capitaine, N., Chapront, J., Lambert, S. and Wallace, P.,
**     "Expressions for the Celestial Intermediate Pole and Celestial
**     Ephemeris Origin consistent with the IAU 2000A precession-
**     nutation model", Astron.Astrophys. 400, 1145-1154 (2003)
**
**     n.b. The celestial ephemeris origin (CEO) was renamed "celestial
**          intermediate origin" (CIO) by IAU 2006 Resolution 2.
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauS06a(date1,date2)
    return s

#iauSp00
try:
    _sofa.iauSp00.argtypes = [
        c_double, #date1
        c_double, #date2
        ]
    _sofa.iauSp00.restype = c_double
except  AttributeError:
    pass
def sp00(date1,date2):
    """
    The TIO locator s', positioning the Terrestrial Intermediate Origin  on the equator of the Celestial Intermediate Pole.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    s : float
        the TIO locator s' in radians (Note 2).
/*
**  - - - - - - - -
**   i a u S p 0 0
**  - - - - - - - -
**
**  The TIO locator s', positioning the Terrestrial Intermediate Origin
**  on the equator of the Celestial Intermediate Pole.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double    TT as a 2-part Julian Date (Note 1)
**
**  Returned (function value):
**                  double    the TIO locator s' in radians (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The TIO locator s' is obtained from polar motion observations by
**     numerical integration, and so is in essence unpredictable.
**     However, it is dominated by a secular drift of about
**     47 microarcseconds per century, which is the approximation
**     evaluated by the present function.
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = _sofa.iauSp00(date1,date2)
    return s


#iauXy06
try:
    _sofa.iauXy06.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#x
        POINTER(c_double),#y
        ]
except  AttributeError:
    pass
def xy06(date1,date2):
    """
    

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    x : float
        CIP X,Y coordinates (Note 2).
    y : float
        CIP X,Y coordinates (Note 2).
/*
**  - - - - - - - -
**   i a u X y 0 6
**  - - - - - - - -
**
**  X,Y coordinates of celestial intermediate pole from series based
**  on IAU 2006 precession and IAU 2000A nutation.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     date1,date2  double     TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double     CIP X,Y coordinates (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The X,Y coordinates are those of the unit vector towards the
**     celestial intermediate pole.  They represent the combined effects
**     of frame bias, precession and nutation.
**
**  3) The fundamental arguments used are as adopted in IERS Conventions
**     (2003) and are from Simon et al. (1994) and Souchay et al.
**     (1999).
**
**  4) This is an alternative to the angles-based method, via the SOFA
**     function iauFw2xy and as used in iauXys06a for example.  The two
**     methods agree at the 1 microarcsecond level (at present), a
**     negligible amount compared with the intrinsic accuracy of the
**     models.  However, it would be unwise to mix the two methods
**     (angles-based and series-based) in a single application.
**
**  Called:
**     iauFal03     mean anomaly of the Moon
**     iauFalp03    mean anomaly of the Sun
**     iauFaf03     mean argument of the latitude of the Moon
**     iauFad03     mean elongation of the Moon from the Sun
**     iauFaom03    mean longitude of the Moon's ascending node
**     iauFame03    mean longitude of Mercury
**     iauFave03    mean longitude of Venus
**     iauFae03     mean longitude of Earth
**     iauFama03    mean longitude of Mars
**     iauFaju03    mean longitude of Jupiter
**     iauFasa03    mean longitude of Saturn
**     iauFaur03    mean longitude of Uranus
**     iauFane03    mean longitude of Neptune
**     iauFapa03    general accumulated precession in longitude
**
**  References:
**
**     Capitaine, N., Wallace, P.T. & Chapront, J., 2003,
**     Astron.Astrophys., 412, 567
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     McCarthy, D. D., Petit, G. (eds.), 2004, IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG
**
**     Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M., 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2019 June 23
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    _sofa.iauXy06(date1,date2,byref(x),byref(y))
    return x.value, y.value


#iauXys00a
try:
    _sofa.iauXys00a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#x
        POINTER(c_double),#y
        POINTER(c_double),#s
        ]
except  AttributeError:
    pass
def xys00a(date1,date2):
    """
    For a given TT date, compute the X,Y coordinates of the Celestial  Intermediate Pole and the CIO locator s, using the IAU 2000A precession-nutation model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).
    s : float
        the CIO locator s (Note 2).
/*
**  - - - - - - - - - -
**   i a u X y s 0 0 a
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2000A
**  precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double   Celestial Intermediate Pole (Note 2)
**     s            double   the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) A faster, but slightly less accurate result (about 1 mas for
**     X,Y), can be obtained by using instead the iauXys00b function.
**
**  Called:
**     iauPnm00a    classical NPB matrix, IAU 2000A
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    s = c_double()
    _sofa.iauXys00a(date1,date2,byref(x),byref(y),byref(s))
    return x.value, y.value, s.value


#iauXys00b
try:
    _sofa.iauXys00b.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#x
        POINTER(c_double),#y
        POINTER(c_double),#s
        ]
except  AttributeError:
    pass
def xys00b(date1,date2):
    """
  For a given TT date, compute the X,Y coordinates of the Celestial  Intermediate Pole and the CIO locator s, using the IAU 2000B precession-nutation model.  

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).
    s : float
        the CIO locator s (Note 2).
/*
**  - - - - - - - - - -
**   i a u X y s 0 0 b
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2000B
**  precession-nutation model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double   TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double   Celestial Intermediate Pole (Note 2)
**     s            double   the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y
**     components of the unit vector in the Geocentric Celestial
**     Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) The present function is faster, but slightly less accurate (about
**     1 mas in X,Y), than the iauXys00a function.
**
**  Called:
**     iauPnm00b    classical NPB matrix, IAU 2000B
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS00       the CIO locator s, given X,Y, IAU 2000A
**
**  Reference:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    s = c_double()
    _sofa.iauXys00b(date1,date2,byref(x),byref(y),byref(s))
    return x.value, y.value, s.value


#iauXys06a
try:
    _sofa.iauXys06a.argtypes = [
        c_double, #date1
        c_double, #date2
        POINTER(c_double),#x
        POINTER(c_double),#y
        POINTER(c_double),#s
        ]
except  AttributeError:
    pass
def xys06a(date1,date2):
    """
    For a given TT date, compute the X,Y coordinates of the Celestial  Intermediate Pole and the CIO locator s, using the IAU 2006 precession and IAU 2000A nutation models.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian Date (Note 1).
    date2 : float
        TT as a 2-part Julian Date (Note 1).

    Returns
    -------
    x : float
        Celestial Intermediate Pole (Note 2).
    y : float
        Celestial Intermediate Pole (Note 2).
    s : float
        the CIO locator s (Note 2).
/*
**  - - - - - - - - - -
**   i a u X y s 0 6 a
**  - - - - - - - - - -
**
**  For a given TT date, compute the X,Y coordinates of the Celestial
**  Intermediate Pole and the CIO locator s, using the IAU 2006
**  precession and IAU 2000A nutation models.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double  TT as a 2-part Julian Date (Note 1)
**
**  Returned:
**     x,y          double  Celestial Intermediate Pole (Note 2)
**     s            double  the CIO locator s (Note 2)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The Celestial Intermediate Pole coordinates are the x,y components
**     of the unit vector in the Geocentric Celestial Reference System.
**
**  3) The CIO locator s (in radians) positions the Celestial
**     Intermediate Origin on the equator of the CIP.
**
**  4) Series-based solutions for generating X and Y are also available:
**     see Capitaine & Wallace (2006) and iauXy06.
**
**  Called:
**     iauPnm06a    classical NPB matrix, IAU 2006/2000A
**     iauBpn2xy    extract CIP X,Y coordinates from NPB matrix
**     iauS06       the CIO locator s, given X,Y, IAU 2006
**
**  References:
**
**     Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
**
**     Wallace, P.T. & Capitaine, N., 2006, Astron.Astrophys. 459, 981
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    x = c_double()
    y = c_double()
    s = c_double()
    _sofa.iauXys06a(date1,date2,byref(x),byref(y),byref(s))
    return x.value, y.value, s.value



#########################################################################
#####   Fundamental arguments for nutation etc.
##########################################################################

#iauFad03
try:
    _sofa.iauFad03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFad03.restype = c_double #a
except  AttributeError:
    pass
def fad03(t):
    """
    
Fundamental argument, IERS Conventions (2003):  mean elongation of the Moon from the Sun.
    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        TDB, Julian centuries since J2000.0 (Note 1).
/*
**  - - - - - - - - -
**   i a u F a d 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean elongation of the Moon from the Sun.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    D, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFad03(t)
    return a


#iauFae03
try:
    _sofa.iauFae03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFae03.restype = c_double #a
except  AttributeError:
    pass
def fae03(t):
    """
     Fundamental argument, IERS Conventions (2003): mean longitude of Earth.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Earth, radians (Note 2).
/*
**  - - - - - - - - -
**   i a u F a e 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Earth.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Earth, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFae03(t)
    return a


#iauFaf03
try:
    _sofa.iauFaf03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFaf03.restype = c_double #a
except  AttributeError:
    pass
def faf03(t):
    """
  Fundamental argument, IERS Conventions (2003):  mean longitude of the Moon minus mean longitude of the ascending  node.  

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        F, radians (Note 2).
/*
**  - - - - - - - - -
**   i a u F a f 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon minus mean longitude of the ascending
**  node.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    F, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFaf03(t)
    return a


#iauFaju03
try:
    _sofa.iauFaju03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFaju03.restype = c_double #a
except  AttributeError:
    pass
def faju03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean longitude of Jupiter.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Jupiter, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a j u 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Jupiter.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Jupiter, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFaju03(t)
    return a


#iauFal03
try:
    _sofa.iauFal03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFal03.restype = c_double #a
except  AttributeError:
    pass
def fal03(t):
    """
     Fundamental argument, IERS Conventions (2003):  mean anomaly of the Moon.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        l, radians (Note 2).
/*
**  - - - - - - - - -
**   i a u F a l 0 3
**  - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Moon.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFal03(t)
    return a


#iauFalp03
try:
    _sofa.iauFalp03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFalp03.restype = c_double #a
except  AttributeError:
    pass
def falp03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean anomaly of the Sun.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        l', radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a l p 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean anomaly of the Sun.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    l', radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFalp03(t)
    return a


#iauFama03
try:
    _sofa.iauFama03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFama03.restype = c_double #a
except  AttributeError:
    pass
def fama03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean longitude of Mars.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Mars, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a m a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mars.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mars, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFama03(t)
    return a


#iauFame03
try:
    _sofa.iauFame03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFame03.restype = c_double #a
except  AttributeError:
    pass
def fame03(t):
    """
    Fundamental argument, IERS Conventions (2003):  mean longitude of Mercury.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Mercury, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a m e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Mercury.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Mercury, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFame03(t)
    return a


#iauFane03
try:
    _sofa.iauFane03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFane03.restype = c_double #a
except  AttributeError:
    pass
def fane03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean longitude of Neptune.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Neptune, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a n e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Neptune.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Neptune, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is adapted from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFane03(t)
    return a


#iauFaom03
try:
    _sofa.iauFaom03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFaom03.restype = c_double #a
except  AttributeError:
    pass
def faom03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean longitude of the Moon's ascending node.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        Omega, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a o m 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of the Moon's ascending node.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    Omega, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFaom03(t)
    return a


#iauFapa03
try:
    _sofa.iauFapa03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFapa03.restype = c_double #a
except  AttributeError:
    pass
def fapa03(t):
    """
    Fundamental argument, IERS Conventions (2003): general accumulated precession in longitude.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        general precession in longitude, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a p a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  general accumulated precession in longitude.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    general precession in longitude, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003).  It
**     is taken from Kinoshita & Souchay (1990) and comes originally
**     from Lieske et al. (1977).
**
**  References:
**
**     Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
**     48, 187
**
**     Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
**     Astron.Astrophys. 58, 1-16
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFapa03(t)
    return a


#iauFasa03
try:
    _sofa.iauFasa03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFasa03.restype = c_double #a
except  AttributeError:
    pass
def fasa03(t):
    """
    Fundamental argument, IERS Conventions (2003):  mean longitude of Saturn.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Saturn, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a s a 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Saturn.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Saturn, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFasa03(t)
    return a


#iauFaur03
try:
    _sofa.iauFaur03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFaur03.restype = c_double #a
except  AttributeError:
    pass
def faur03(t):
    """
    
Fundamental argument, IERS Conventions (2003):  mean longitude of Uranus.
    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Uranus, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a u r 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Uranus.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned  (function value):
**           double    mean longitude of Uranus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     is adapted from Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFaur03(t)
    return a


#iauFave03
try:
    _sofa.iauFave03.argtypes = [
        c_double, #t
        ]
    _sofa.iauFave03.restype = c_double #a
except  AttributeError:
    pass
def fave03(t):
    """
    Fundamental argument, IERS Conventions (2003): mean longitude of Venus.

    Parameters
    ----------
    t : float
        TDB, Julian centuries since J2000.0 (Note 1).

    Returns
    -------
    a : float
        mean longitude of Venus, radians (Note 2).
/*
**  - - - - - - - - - -
**   i a u F a v e 0 3
**  - - - - - - - - - -
**
**  Fundamental argument, IERS Conventions (2003):
**  mean longitude of Venus.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  canonical model.
**
**  Given:
**     t     double    TDB, Julian centuries since J2000.0 (Note 1)
**
**  Returned (function value):
**           double    mean longitude of Venus, radians (Note 2)
**
**  Notes:
**
**  1) Though t is strictly TDB, it is usually more convenient to use
**     TT, which makes no significant difference.
**
**  2) The expression used is as adopted in IERS Conventions (2003) and
**     comes from Souchay et al. (1999) after Simon et al. (1994).
**
**  References:
**
**     McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
**     IERS Technical Note No. 32, BKG (2004)
**
**     Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
**     Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683
**
**     Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
**     Astron.Astrophys.Supp.Ser. 135, 111
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = _sofa.iauFave03(t)
    return a


#############################################################################
####  Star catalog conversions
#############################################################################

#iauFk45z
try:
    _sofa.iauFk45z.argtypes = [
        c_double, # r1950
        c_double, # d1950
        c_double, # bepoch
        POINTER(c_double), # r2000
        POINTER(c_double), # d2000
        ]
except  AttributeError:
    pass
def fk45z(r1950,d1950,bepoch):
    """
    Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero proper motion in the FK5 system.

    Parameters
    ----------
    r1950 : float
        B1950.0 FK4 RA,Dec at epoch (rad).
    d1950 : float
        B1950.0 FK4 RA,Dec at epoch (rad).
    bepoch : float
        Besselian epoch (e.g. 1979.3D0).

    Returns
    -------
    r2000 : float
        J2000.0 FK5 RA,Dec (rad).
    d2000 : float
        J2000.0 FK5 RA,Dec (rad).
/*
**  - - - - - - - - -
**   i a u F k 4 5 z
**  - - - - - - - - -
**
**  Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
**  proper motion in the FK5 system.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  This function converts a star's catalog data from the old FK4
**  (Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system,
**  in such a way that the FK5 proper motion is zero.  Because such a
**  star has, in general, a non-zero proper motion in the FK4 system,
**  the routine requires the epoch at which the position in the FK4
**  system was determined.
**
**  Given:
**     r1950,d1950    double   B1950.0 FK4 RA,Dec at epoch (rad)
**     bepoch         double   Besselian epoch (e.g. 1979.3D0)
**
**  Returned:
**     r2000,d2000    double   J2000.0 FK5 RA,Dec (rad)
**
**  Notes:
**
**  1) The epoch bepoch is strictly speaking Besselian, but if a
**     Julian epoch is supplied the result will be affected only to a
**     negligible extent.
**
**  2) The method is from Appendix 2 of Aoki et al. (1983), but using
**     the constants of Seidelmann (1992).  See the routine iauFk425
**     for a general introduction to the FK4 to FK5 conversion.
**
**  3) Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only
**     is provided for.  Conversions for different starting and/or
**     ending epochs would require additional treatment for precession,
**     proper motion and E-terms.
**
**  4) In the FK4 catalog the proper motions of stars within 10 degrees
**     of the poles do not embody differential E-terms effects and
**     should, strictly speaking, be handled in a different manner from
**     stars outside these regions.  However, given the general lack of
**     homogeneity of the star data available for routine astrometry,
**     the difficulties of handling positions that may have been
**     determined from astrometric fields spanning the polar and non-
**     polar regions, the likelihood that the differential E-terms
**     effect was not taken into account when allowing for proper motion
**     in past astrometry, and the undesirability of a discontinuity in
**     the algorithm, the decision has been made in this SOFA algorithm
**     to include the effects of differential E-terms on the proper
**     motions for all stars, whether polar or not.  At epoch 2000.0,
**     and measuring "on the sky" rather than in terms of RA change, the
**     errors resulting from this simplification are less than
**     1 milliarcsecond in position and 1 milliarcsecond per century in
**     proper motion.
**
**  References:
**
**     Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
**     FK4-based positions of stars to epoch J2000.0 positions in
**     accordance with the new IAU resolutions".  Astron.Astrophys.
**     128, 263-267.
**
**     Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
**     Astronomical Almanac", ISBN 0-935702-68-7.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauC2s       p-vector to spherical
**     iauEpb2jd    Besselian epoch to Julian date
**     iauEpj       Julian date to Julian epoch
**     iauPdp       scalar product of two p-vectors
**     iauPmp       p-vector minus p-vector
**     iauPpsp      p-vector plus scaled p-vector
**     iauPvu       update a pv-vector
**     iauS2c       spherical to p-vector
**
**  This revision:   2018 December 5
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r2000 = c_double()
    d2000 = c_double()
    _sofa.iauFk45z(r1950,d1950,bepoch,byref(r2000),byref(d2000))
    return r2000.value, d2000.value


#iauFk425
try:
    _sofa.iauFk425.argtypes = [
        c_double, # r1950
        c_double, # d1950
        c_double, # dr1950
        c_double, # dd1950
        c_double, # p1950
        c_double, # v1950
        POINTER(c_double), # r2000
        POINTER(c_double), # d2000
        POINTER(c_double), # dr2000
        POINTER(c_double), # dd2000
        POINTER(c_double), # p2000
        POINTER(c_double), # v2000
        ]
except  AttributeError:
    pass
def fk425(r1950,d1950,dr1950,dd1950,p1950,v1950):
    """
    Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

    Parameters
    ----------
    r1950 : flaot
        B1950.0 RA,Dec (rad).
    d1950 : flaot
        B1950.0 RA,Dec (rad).
    dr1950 : flaot
        B1950.0 proper motions (rad/trop.yr).
    dd1950 : flaot
        B1950.0 proper motions (rad/trop.yr).
    p1950 : flaot
        parallax (arcsec).
    v1950 : flaot
        radial velocity (km/s, +ve = moving away).

    Returns
    -------
    r2000 : flaot
        J2000.0 RA,Dec (rad).
    d2000 : flaot
        J2000.0 RA,Dec (rad).
    dr2000 : flaot
        J2000.0 proper motions (rad/Jul.yr).
    dd2000 : flaot
        J2000.0 proper motions (rad/Jul.yr).
    p2000 : flaot
        parallax (arcsec).
    v2000 : TYPE
        radial velocity (km/s, +ve = moving away).
/*
**  - - - - - - - - -
**   i a u F k 4 2 5
**  - - - - - - - - -
**
**  Convert B1950.0 FK4 star catalog data to J2000.0 FK5.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  This function converts a star's catalog data from the old FK4
** (Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system.
**
**  Given: (all B1950.0, FK4)
**     r1950,d1950    double   B1950.0 RA,Dec (rad)
**     dr1950,dd1950  double   B1950.0 proper motions (rad/trop.yr)
**     p1950          double   parallax (arcsec)
**     v1950          double   radial velocity (km/s, +ve = moving away)
**
**  Returned: (all J2000.0, FK5)
**     r2000,d2000    double   J2000.0 RA,Dec (rad)
**     dr2000,dd2000  double   J2000.0 proper motions (rad/Jul.yr)
**     p2000          double   parallax (arcsec)
**     v2000          double   radial velocity (km/s, +ve = moving away)
**
**  Notes:
**
**  1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
**     and are per year rather than per century.
**
**  2) The conversion is somewhat complicated, for several reasons:
**
**     . Change of standard epoch from B1950.0 to J2000.0.
**
**     . An intermediate transition date of 1984 January 1.0 TT.
**
**     . A change of precession model.
**
**     . Change of time unit for proper motion (tropical to Julian).
**
**     . FK4 positions include the E-terms of aberration, to simplify
**       the hand computation of annual aberration.  FK5 positions
**       assume a rigorous aberration computation based on the Earth's
**       barycentric velocity.
**
**     . The E-terms also affect proper motions, and in particular cause
**       objects at large distances to exhibit fictitious proper
**       motions.
**
**     The algorithm is based on Smith et al. (1989) and Yallop et al.
**     (1989), which presented a matrix method due to Standish (1982) as
**     developed by Aoki et al. (1983), using Kinoshita's development of
**     Andoyer's post-Newcomb precession.  The numerical constants from
**     Seidelmann (1992) are used canonically.
**
**  3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
**     Conversions for different epochs and equinoxes would require
**     additional treatment for precession, proper motion and E-terms.
**
**  4) In the FK4 catalog the proper motions of stars within 10 degrees
**     of the poles do not embody differential E-terms effects and
**     should, strictly speaking, be handled in a different manner from
**     stars outside these regions.  However, given the general lack of
**     homogeneity of the star data available for routine astrometry,
**     the difficulties of handling positions that may have been
**     determined from astrometric fields spanning the polar and non-
**     polar regions, the likelihood that the differential E-terms
**     effect was not taken into account when allowing for proper motion
**     in past astrometry, and the undesirability of a discontinuity in
**     the algorithm, the decision has been made in this SOFA algorithm
**     to include the effects of differential E-terms on the proper
**     motions for all stars, whether polar or not.  At epoch J2000.0,
**     and measuring "on the sky" rather than in terms of RA change, the
**     errors resulting from this simplification are less than
**     1 milliarcsecond in position and 1 milliarcsecond per century in
**     proper motion.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauPv2s      pv-vector to spherical coordinates
**     iauPdp       scalar product of two p-vectors
**     iauPvmpv     pv-vector minus pv_vector
**     iauPvppv     pv-vector plus pv_vector
**     iauS2pv      spherical coordinates to pv-vector
**     iauSxp       multiply p-vector by scalar
**
**  References:
**
**     Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
**     FK4-based positions of stars to epoch J2000.0 positions in
**     accordance with the new IAU resolutions".  Astron.Astrophys.
**     128, 263-267.
**
**     Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
**     Astronomical Almanac", ISBN 0-935702-68-7.
**
**     Smith, C.A. et al., 1989, "The transformation of astrometric
**     catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
**
**     Standish, E.M., 1982, "Conversion of positions and proper motions
**     from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
**     115, 1, 20-22.
**
**     Yallop, B.D. et al., 1989, "Transformation of mean star places
**     from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
**     Astron.J. 97, 274.
**
**  This revision:   2018 December 5
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r2000 = c_double()
    d2000 = c_double()
    dr2000 = c_double()
    dd2000 = c_double()
    p2000 = c_double()
    v2000 = c_double()
    _sofa.iauFk425(r1950,d1950,dr1950,dd1950,p1950,v1950,byref(r2000),byref(d2000),byref(dr2000),byref(dd2000),byref(p2000),byref(v2000))
    return r2000.value, d2000.value, dr2000.value, dd2000.value,p2000.value,v2000.value


#iauFk5hip
try:
    _sofa.iauFk5hip.argtypes = [
        ndpointer(shape = (3,3),dtype = float),#r5h
        ndpointer(shape = 3,dtype = float),#s5h
        ]
except  AttributeError:
    pass
def fk5hip():
    """
    FK5 to Hipparcos rotation and spin.

    Returns
    -------
    r5h : double[3][3]
        r-matrix: FK5 rotation wrt Hipparcos (Note 2).
    s5h : double[3][3]
        r-vector: FK5 spin wrt Hipparcos (Note 3).
/*
**  - - - - - - - - - -
**   i a u F k 5 h i p
**  - - - - - - - - - -
**
**  FK5 to Hipparcos rotation and spin.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Returned:
**     r5h   double[3][3]  r-matrix: FK5 rotation wrt Hipparcos (Note 2)
**     s5h   double[3]     r-vector: FK5 spin wrt Hipparcos (Note 3)
**
**  Notes:
**
**  1) This function models the FK5 to Hipparcos transformation as a
**     pure rotation and spin;  zonal errors in the FK5 catalogue are
**     not taken into account.
**
**  2) The r-matrix r5h operates in the sense:
**
**           P_Hipparcos = r5h x P_FK5
**
**     where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
**     the equivalent Hipparcos p-vector.
**
**  3) The r-vector s5h represents the time derivative of the FK5 to
**     Hipparcos rotation.  The units are radians per year (Julian,
**     TDB).
**
**  Called:
**     iauRv2m      r-vector to r-matrix
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
**
**  This revision:  2017 October 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r5h = zeros(shape = (3,3),dtype = float)
    s5h = zeros(shape = 3,dtype = float)
    _sofa.iauFk5hip(r5h,s5h)
    return r5h,s5h


#iauFk5hz
try:
    _sofa.iauFk5hz.argtypes = [
        c_double, # r5
        c_double, # d5
        c_double, # date1
        c_double, # date2
        POINTER(c_double), # rh
        POINTER(c_double), # dh
        ]
except  AttributeError:
    pass
def fk5hz(r5,d5,date1,date2):
    """
    Transform an FK5 (J2000.0) star position into the system of the Hipparcos catalogue, assuming zero Hipparcos proper motion.

    Parameters
    ----------
    r5 : float
        FK5 RA (radians), equinox J2000.0, at date.
    d5 : float
         FK5 Dec (radians), equinox J2000.0, at date.
    date1 : float
        TDB date (Notes 1,2).
    date2 : float
        TDB date (Notes 1,2).

    Returns
    -------
   rh : float
        Hipparcos RA (radians).
    dh : float
        Hipparcos Dec (radians).
/*
**  - - - - - - - - -
**   i a u F k 5 h z
**  - - - - - - - - -
**
**  Transform an FK5 (J2000.0) star position into the system of the
**  Hipparcos catalogue, assuming zero Hipparcos proper motion.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     r5           double   FK5 RA (radians), equinox J2000.0, at date
**     d5           double   FK5 Dec (radians), equinox J2000.0, at date
**     date1,date2  double   TDB date (Notes 1,2)
**
**  Returned:
**     rh           double   Hipparcos RA (radians)
**     dh           double   Hipparcos Dec (radians)
**
**  Notes:
**
**  1) This function converts a star position from the FK5 system to
**     the Hipparcos system, in such a way that the Hipparcos proper
**     motion is zero.  Because such a star has, in general, a non-zero
**     proper motion in the FK5 system, the function requires the date
**     at which the position in the FK5 system was determined.
**
**  2) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalogue are not
**     taken into account.
**
**  4) The position returned by this function is in the Hipparcos
**     reference system but at date date1+date2.
**
**  5) See also iauFk52h, iauH2fk5, iauHfk5z.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauFk5hip    FK5 to Hipparcos rotation and spin
**     iauSxp       multiply p-vector by scalar
**     iauRv2m      r-vector to r-matrix
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauPxp       vector product of two p-vectors
**     iauC2s       p-vector to spherical
**     iauAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rh = c_double()
    dh = c_double()
    _sofa.iauFk5hz(r5,d5,date1,date2,byref(rh),byref(dh))
    return rh.value, dh.value


#iauFk52h
try:
    _sofa.iauFk52h.argtypes = [
        c_double, # r5
        c_double, # d5
        c_double, # dr5
        c_double, # dd5
        c_double, # px5
        c_double, # rv5
        POINTER(c_double), # rh
        POINTER(c_double), # dh
        POINTER(c_double), # drh
        POINTER(c_double), # ddh
        POINTER(c_double), # pxh
        POINTER(c_double), # rvh
        ]
except AttributeError:
    pass
def fk52h(r5,d5,dr5,dd5,px5,rv5):
    """
    Transform FK5 (J2000.0) star data into the Hipparcos system.

    Parameters
    ----------
    r5 : float
        RA (radians).
    d5 : float
        Dec (radians).
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear).
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear).
    px5 : float
        parallax (arcsec).
    rv5 : float
        radial velocity (km/s, positive = receding).

    Returns
    -------
    rh : float
        RA (radians).
    dh : float
        Dec (radians).
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear).
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear).
    pxh : float
        parallax (arcsec).
    rvh : float
        radial velocity (km/s, positive = receding).
/*
**  - - - - - - - - -
**   i a u F k 5 2 h
**  - - - - - - - - -
**
**  Transform FK5 (J2000.0) star data into the Hipparcos system.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given (all FK5, equinox J2000.0, epoch J2000.0):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     px5     double    parallax (arcsec)
**     rv5     double    radial velocity (km/s, positive = receding)
**
**  Returned (all Hipparcos, epoch J2000.0):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     pxh     double    parallax (arcsec)
**     rvh     double    radial velocity (km/s, positive = receding)
**
**  Notes:
**
**  1) This function transforms FK5 star positions and proper motions
**     into the system of the Hipparcos catalog.
**
**  2) The proper motions in RA are dRA/dt rather than
**     cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalog are not
**     taken into account.
**
**  4) See also iauH2fk5, iauFk5hz, iauHfk5z.
**
**  Called:
**     iauStarpv    star catalog data to space motion pv-vector
**     iauFk5hip    FK5 to Hipparcos rotation and spin
**     iauRxp       product of r-matrix and p-vector
**     iauPxp       vector product of two p-vectors
**     iauPpp       p-vector plus p-vector
**     iauPvstar    space motion pv-vector to star catalog data
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
**
**  This revision:  2017 October 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rh = c_double()
    dh = c_double()
    drh = c_double()
    ddh = c_double()
    pxh = c_double()
    rvh = c_double()
    _sofa.iauFk52h(r5,d5,dr5,dd5,px5,rv5,byref(rh),byref(dh),byref(drh),byref(ddh),byref(pxh),byref(rvh))
    return rh.value, dh.value, drh.value, ddh.value, pxh.value, rvh.value


#iauFk54z
try:
    _sofa.iauFk54z.argtypes = [
        c_double, # r2000
        c_double, # d2000
        c_double, # bepoch
        POINTER(c_double), # r1950
        POINTER(c_double), # d1950
        POINTER(c_double), # dr1950
        POINTER(c_double), # dd1950
        ]
except AttributeError:
    pass
def fk54z(r2000,d2000,bepoch):
    """
    Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero  proper motion in FK5 and parallax.v

    Parameters
    ----------
    r2000 : float
        J2000.0 FK5 RA,Dec (rad).
    d2000 : float
        J2000.0 FK5 RA,Dec (rad).
    bepoch : float
        Besselian epoch (e.g. 1950.0).

    Returns
    -------
    r1950 : float
        B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH.
    d1950 : float
        B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH.
    dr1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr).
    dd1950 : float
        B1950.0 FK4 proper motions (rad/trop.yr).
/*
**  - - - - - - - - -
**   i a u F k 5 4 z
**  - - - - - - - - -
**
**  Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
**  proper motion in FK5 and parallax.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     r2000,d2000    double   J2000.0 FK5 RA,Dec (rad)
**     bepoch         double   Besselian epoch (e.g. 1950.0)
**
**  Returned:
**     r1950,d1950    double   B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH
**     dr1950,dd1950  double   B1950.0 FK4 proper motions (rad/trop.yr)
**
**  Notes:
**
**  1) In contrast to the iauFk524  routine, here the FK5 proper
**     motions, the parallax and the radial velocity are presumed zero.
**
**  2) This function converts a star position from the IAU 1976 FK5
**    (Fricke) system to the former FK4 (Bessel-Newcomb) system, for
**     cases such as distant radio sources where it is presumed there is
**     zero parallax and no proper motion.  Because of the E-terms of
**     aberration, such objects have (in general) non-zero proper motion
**     in FK4, and the present routine returns those fictitious proper
**     motions.
**
**  3) Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
**     Conversions involving other equinoxes would require additional
**     treatment for precession.
**
**  4) The position returned by this routine is in the B1950.0 FK4
**     reference system but at Besselian epoch BEPOCH.  For comparison
**     with catalogs the BEPOCH argument will frequently be 1950.0. (In
**     this context the distinction between Besselian and Julian epoch
**     is insignificant.)
**
**  5) The RA component of the returned (fictitious) proper motion is
**     dRA/dt rather than cos(Dec)*dRA/dt.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauC2s       p-vector to spherical
**     iauFk524     FK4 to FK5
**     iauS2c       spherical to p-vector
**
**  This revision:   2018 December 5
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r1950 = c_double()
    d1950 = c_double()
    dr1950 = c_double()
    dd1950 = c_double()
    _sofa.iauFk54z(r2000,d2000,bepoch,byref(r1950),byref(d1950),byref(dr1950),byref(dd1950))
    return r1950.value, d1950.value, dr1950.value, dd1950.value


#iauFk524
try:
    _sofa.iauFk524.argtypes = [
        c_double, # r2000
        c_double, # d2000
        c_double, # dr2000
        c_double, # dd2000
        c_double, #p2000
        c_double, # v2000
        POINTER(c_double), # r1950
        POINTER(c_double), # d1950
        POINTER(c_double), # dr1950
        POINTER(c_double), # dd1950
        POINTER(c_double), # p1950
        POINTER(c_double), # v1950
        ]
except AttributeError:
    pass
def fk524(r2000,d2000,dr2000,dd2000,p2000,v2000):
    """
    Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

    Parameters
    ----------
    r2000 : float
        J2000.0 RA,Dec (rad).
    d2000 : float
        J2000.0 RA,Dec (rad).
    dr2000 : float
        J2000.0 proper motions (rad/Jul.yr).
    dd2000 : float
        J2000.0 proper motions (rad/Jul.yr).
    p2000 : float
        parallax (arcsec).
    v2000 : float
        radial velocity (km/s, +ve = moving away).

    Returns
    -------
    r1950 : float
        B1950.0 RA,Dec (rad).
    d1950 : float
        B1950.0 RA,Dec (rad).
    dr1950 : float
        B1950.0 proper motions (rad/trop.yr).
    dd1950 : float
        B1950.0 proper motions (rad/trop.yr).
    p1950 : float
        parallax (arcsec).
    v1950 : float
        radial velocity (km/s, +ve = moving away).
/*
**  - - - - - - - - -
**   i a u F k 5 2 4
**  - - - - - - - - -
**
**  Convert J2000.0 FK5 star catalog data to B1950.0 FK4.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given: (all J2000.0, FK5)
**     r2000,d2000    double   J2000.0 RA,Dec (rad)
**     dr2000,dd2000  double   J2000.0 proper motions (rad/Jul.yr)
**     p2000          double   parallax (arcsec)
**     v2000          double   radial velocity (km/s, +ve = moving away)
**
**  Returned: (all B1950.0, FK4)
**     r1950,d1950    double   B1950.0 RA,Dec (rad)
**     dr1950,dd1950  double   B1950.0 proper motions (rad/trop.yr)
**     p1950          double   parallax (arcsec)
**     v1950          double   radial velocity (km/s, +ve = moving away)
**
**  Notes:
**
**  1) The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
**     and are per year rather than per century.
**
**  2) The conversion is somewhat complicated, for several reasons:
**
**     . Change of standard epoch from J2000.0 to B1950.0.
**
**     . An intermediate transition date of 1984 January 1.0 TT.
**
**     . A change of precession model.
**
**     . Change of time unit for proper motion (Julian to tropical).
**
**     . FK4 positions include the E-terms of aberration, to simplify
**       the hand computation of annual aberration.  FK5 positions
**       assume a rigorous aberration computation based on the Earth's
**       barycentric velocity.
**
**     . The E-terms also affect proper motions, and in particular cause
**       objects at large distances to exhibit fictitious proper
**       motions.
**
**     The algorithm is based on Smith et al. (1989) and Yallop et al.
**     (1989), which presented a matrix method due to Standish (1982) as
**     developed by Aoki et al. (1983), using Kinoshita's development of
**     Andoyer's post-Newcomb precession.  The numerical constants from
**     Seidelmann (1992) are used canonically.
**
**  4) In the FK4 catalog the proper motions of stars within 10 degrees
**     of the poles do not embody differential E-terms effects and
**     should, strictly speaking, be handled in a different manner from
**     stars outside these regions.  However, given the general lack of
**     homogeneity of the star data available for routine astrometry,
**     the difficulties of handling positions that may have been
**     determined from astrometric fields spanning the polar and non-
**     polar regions, the likelihood that the differential E-terms
**     effect was not taken into account when allowing for proper motion
**     in past astrometry, and the undesirability of a discontinuity in
**     the algorithm, the decision has been made in this SOFA algorithm
**     to include the effects of differential E-terms on the proper
**     motions for all stars, whether polar or not.  At epoch J2000.0,
**     and measuring "on the sky" rather than in terms of RA change, the
**     errors resulting from this simplification are less than
**     1 milliarcsecond in position and 1 milliarcsecond per century in
**     proper motion.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauPdp       scalar product of two p-vectors
**     iauPm        modulus of p-vector
**     iauPmp       p-vector minus p-vector
**     iauPpp       p-vector pluus p-vector
**     iauPv2s      pv-vector to spherical coordinates
**     iauS2pv      spherical coordinates to pv-vector
**     iauSxp       multiply p-vector by scalar
**
**  References:
**
**     Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
**     FK4-based positions of stars to epoch J2000.0 positions in
**     accordance with the new IAU resolutions".  Astron.Astrophys.
**     128, 263-267.
**
**     Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
**     Astronomical Almanac", ISBN 0-935702-68-7.
**
**     Smith, C.A. et al., 1989, "The transformation of astrometric
**     catalog systems to the equinox J2000.0".  Astron.J. 97, 265.
**
**     Standish, E.M., 1982, "Conversion of positions and proper motions
**     from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
**     115, 1, 20-22.
**
**     Yallop, B.D. et al., 1989, "Transformation of mean star places
**     from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
**     Astron.J. 97, 274.
**
**  This revision:   2018 December 5
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r1950 = c_double()
    d1950 = c_double()
    dr1950 = c_double()
    dd1950= c_double()
    p1950 = c_double()
    v1950 = c_double()
    _sofa.iauFk524(r2000,d2000,dr2000,dd2000,p2000,v2000,byref(r1950),byref(d1950),byref(dr1950),byref(dd1950),byref(p1950),byref(v1950))
    return r1950.value, d1950.value, dr1950.value, dd1950.value, p1950.value, v1950.value


#iauH2fk5
try:
    _sofa.iauH2fk5.argtypes = [
        c_double, # rh
        c_double, # dh
        c_double, # drh
        c_double, # ddh
        c_double, # pxh
        c_double, #rvh
        POINTER(c_double), # r5
        POINTER(c_double), # d5
        POINTER(c_double), # dr5
        POINTER(c_double), # dd5
        POINTER(c_double), # px5
        POINTER(c_double), # rv5
        ]
except AttributeError:
    pass
def h2fk5(rh,dh,drh,ddh,pxh,rvh):
    """
    Transform Hipparcos star data into the FK5 (J2000.0) system.

    Parameters
    (all Hipparcos, epoch J2000.0)
    ----------
    rh : float
        RA (radians).
    dh : float
        Dec (radians).
    drh : float
        proper motion in RA (dRA/dt, rad/Jyear).
    ddh : float
        proper motion in Dec (dDec/dt, rad/Jyear).
    pxh : float
        parallax (arcsec).
    rvh : float
        radial velocity (km/s, positive = receding).

    Returns
    (all FK5, equinox J2000.0, epoch J2000.0)
    -------
    r5 : float
        RA (radians).
    d5 : float
        Dec (radians).
    dr5 : float
        proper motion in RA (dRA/dt, rad/Jyear).
    dd5 : float
        proper motion in Dec (dDec/dt, rad/Jyear).
    px5 : float
        parallax (arcsec).
    rv5 : float
        radial velocity (km/s, positive = receding).
/*
**  - - - - - - - - -
**   i a u H 2 f k 5
**  - - - - - - - - -
**
**  Transform Hipparcos star data into the FK5 (J2000.0) system.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given (all Hipparcos, epoch J2000.0):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     pxh     double    parallax (arcsec)
**     rvh     double    radial velocity (km/s, positive = receding)
**
**  Returned (all FK5, equinox J2000.0, epoch J2000.0):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**     px5     double    parallax (arcsec)
**     rv5     double    radial velocity (km/s, positive = receding)
**
**  Notes:
**
**  1) This function transforms Hipparcos star positions and proper
**     motions into FK5 J2000.0.
**
**  2) The proper motions in RA are dRA/dt rather than
**     cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure
**     rotation and spin;  zonal errors in the FK5 catalog are not
**     taken into account.
**
**  4) See also iauFk52h, iauFk5hz, iauHfk5z.
**
**  Called:
**     iauStarpv    star catalog data to space motion pv-vector
**     iauFk5hip    FK5 to Hipparcos rotation and spin
**     iauRv2m      r-vector to r-matrix
**     iauRxp       product of r-matrix and p-vector
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauPxp       vector product of two p-vectors
**     iauPmp       p-vector minus p-vector
**     iauPvstar    space motion pv-vector to star catalog data
**
**  Reference:
**
**     F.Mignard & M.Froeschle, Astron.Astrophys., 354, 732-739 (2000).
**
**  This revision:  2017 October 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r5 = c_double()
    d5 = c_double()
    dr5 = c_double()
    dd5 = c_double()
    px5 = c_double()
    rv5 = c_double()
    _sofa.iauH2fk5(rh,dh,drh,ddh,pxh,rvh,byref(r5),byref(d5),byref(dr5),byref(dd5),byref(px5),byref(rv5))
    return r5.value, d5.value, dr5.value, dd5.value, px5.value, rv5.value


#iauHfk5z
try:
    _sofa.iauHfk5z.argtypes = [
        c_double, # rh
        c_double, # dh
        c_double, # date1
        c_double, # date2
        POINTER(c_double), # r5
        POINTER(c_double), # d5
        POINTER(c_double), # dr5
        POINTER(c_double), # dd5
        ]
except AttributeError:
    pass
def hfk5z(rh,dh,date1,date2):
    """
    Transform a Hipparcos star position into FK5 J2000.0, assuming zero Hipparcos proper motion.

    Parameters
    ----------
    rh : float
        Hipparcos RA (radians).
    dh : float
         Hipparcos Dec (radians).
    date1 : float
        TDB date (Note 1).
    date2 : float
        TDB date (Note 1).

    Returns
    (all FK5, equinox J2000.0, date date1+date2)
    -------
    r5 : float
        RA (radians).
    d5 : float
        Dec (radians).
    dr5 : float
        FK5 RA proper motion (rad/year, Note 4).
    dd5 : float
        Dec proper motion (rad/year, Note 4).
/*
**  - - - - - - - - -
**   i a u H f k 5 z
**  - - - - - - - - -
**
**  Transform a Hipparcos star position into FK5 J2000.0, assuming
**  zero Hipparcos proper motion.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     rh            double    Hipparcos RA (radians)
**     dh            double    Hipparcos Dec (radians)
**     date1,date2   double    TDB date (Note 1)
**
**  Returned (all FK5, equinox J2000.0, date date1+date2):
**     r5            double    RA (radians)
**     d5            double    Dec (radians)
**     dr5           double    FK5 RA proper motion (rad/year, Note 4)
**     dd5           double    Dec proper motion (rad/year, Note 4)
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
**
**  3) The FK5 to Hipparcos transformation is modeled as a pure rotation
**     and spin;  zonal errors in the FK5 catalogue are not taken into
**     account.
**
**  4) It was the intention that Hipparcos should be a close
**     approximation to an inertial frame, so that distant objects have
**     zero proper motion;  such objects have (in general) non-zero
**     proper motion in FK5, and this function returns those fictitious
**     proper motions.
**
**  5) The position returned by this function is in the FK5 J2000.0
**     reference system but at date date1+date2.
**
**  6) See also iauFk52h, iauH2fk5, iauFk5zhz.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauFk5hip    FK5 to Hipparcos rotation and spin
**     iauRxp       product of r-matrix and p-vector
**     iauSxp       multiply p-vector by scalar
**     iauRxr       product of two r-matrices
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauPxp       vector product of two p-vectors
**     iauPv2s      pv-vector to spherical
**     iauAnp       normalize angle into range 0 to 2pi
**
**  Reference:
**
**     F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r5 = c_double()
    d5 = c_double()
    dr5 = c_double()
    dd5 = c_double()
    _sofa.iauHfk5z(rh,dh,date1,date2,byref(r5),byref(d5),byref(dr5),byref(dd5))
    return r5.value, d5.value, dr5.value, dd5.value



############################################################################
#########    Ecliptic coordinates
##########################################################################

#iauEceq06
try:
    _sofa.iauEceq06.argtypes = [
        c_double, # date1
        c_double, # date2
        c_double, # dl
        c_double, # db
        POINTER(c_double), # dr
        POINTER(c_double), # dd
        ]
except AttributeError:
    pass
def eceq06(date1,date2,dl,db):
    """
    Transformation from ecliptic coordinates (mean equinox and ecliptic  of date) to ICRS RA,Dec, using the IAU 2006 precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian date (Note 1).
    date2 : float
        TT as a 2-part Julian date (Note 1).
    dl : float
        ecliptic longitude and latitude (radians).
    db : float
        ecliptic longitude and latitude (radians).

    Returns
    -------
    dr : float
        ICRS right ascension  (radians).
    dd : float
        ICRS declination (radians).
/*
**  - - - - - - - - - -
**   i a u E c e q 0 6
**  - - - - - - - - - -
**
**  Transformation from ecliptic coordinates (mean equinox and ecliptic
**  of date) to ICRS RA,Dec, using the IAU 2006 precession model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double TT as a 2-part Julian date (Note 1)
**     dl,db       double ecliptic longitude and latitude (radians)
**
**  Returned:
**     dr,dd       double ICRS right ascension and declination (radians)
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) No assumptions are made about whether the coordinates represent
**     starlight and embody astrometric effects such as parallax or
**     aberration.
**
**  3) The transformation is approximately that from ecliptic longitude
**     and latitude (mean equinox and ecliptic of date) to mean J2000.0
**     right ascension and declination, with only frame bias (always
**     less than 25 mas) to disturb this classical picture.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauEcm06     J2000.0 to ecliptic rotation matrix, IAU 2006
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauC2s       unit vector to spherical coordinates
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dr = c_double()
    dd = c_double()
    _sofa.iauEceq06(date1,date2,dl,db,byref(dr),byref(dd))
    return dr.value, dd.value


#iauEcm06
try:
    _sofa.iauEcm06.argtypes = [
        c_double, # date1
        c_double, # date2
        ndpointer(shape=(3,3),dtype = float),#rm[3][3]
        ]
except AttributeError:
    pass
def ecm06(date1,date2):
    """
    ICRS equatorial to ecliptic rotation matrix, IAU 2006.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian date (Note 1).
    date2 : float
        TT as a 2-part Julian date (Note 1).

    Returns
    -------
    rm : double[3][3]
        ICRS to ecliptic rotation matrix.
/*
**  - - - - - - - - -
**   i a u E c m 0 6
**  - - - - - - - - -
**
**  ICRS equatorial to ecliptic rotation matrix, IAU 2006.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2  double         TT as a 2-part Julian date (Note 1)
**
**  Returned:
**     rm           double[3][3]   ICRS to ecliptic rotation matrix
**
**  Notes:
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  1) The matrix is in the sense
**
**        E_ep = rm x P_ICRS,
**
**     where P_ICRS is a vector with respect to ICRS right ascension
**     and declination axes and E_ep is the same vector with respect to
**     the (inertial) ecliptic and equinox of date.
**
**  2) P_ICRS is a free vector, merely a direction, typically of unit
**     magnitude, and not bound to any particular spatial origin, such
**     as the Earth, Sun or SSB.  No assumptions are made about whether
**     it represents starlight and embodies astrometric effects such as
**     parallax or aberration.  The transformation is approximately that
**     between mean J2000.0 right ascension and declination and ecliptic
**     longitude and latitude, with only frame bias (always less than
**     25 mas) to disturb this classical picture.
**
**  Called:
**     iauObl06     mean obliquity, IAU 2006
**     iauPmat06    PB matrix, IAU 2006
**     iauIr        initialize r-matrix to identity
**     iauRx        rotate around X-axis
**     iauRxr       product of two r-matrices
**
**  This revision:  2015 December 11
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rm = zeros(shape=(3,3),dtype=float)
    _sofa.iauEcm06(date1,date2,rm)
    return rm


#iauEqec06
try:
    _sofa.iauEqec06.argtypes = [
        c_double, # date1
        c_double, # date2
        c_double, # dr
        c_double, # dd
        POINTER(c_double), # dl
        POINTER(c_double), # db
        ]
except AttributeError:
    pass
def eqec06(date1,date2,dr,dd):
    """
    Transformation from ICRS equatorial coordinates to ecliptic  coordinates (mean equinox and ecliptic of date) using IAU 2006  precession model.

    Parameters
    ----------
    date1 : float
        TT as a 2-part Julian date (Note 1).
    date2 : float
        TT as a 2-part Julian date (Note 1).
    dr : float
        ICRS right ascension and declination (radians).
    dd : float
        ICRS right ascension and declination (radians).

    Returns
    -------
    dl : float
        ecliptic longitude and latitude (radians).
    db : float
        ecliptic longitude and latitude (radians).
/*
**  - - - - - - - - - -
**   i a u E q e c 0 6
**  - - - - - - - - - -
**
**  Transformation from ICRS equatorial coordinates to ecliptic
**  coordinates (mean equinox and ecliptic of date) using IAU 2006
**  precession model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     date1,date2 double TT as a 2-part Julian date (Note 1)
**     dr,dd       double ICRS right ascension and declination (radians)
**
**  Returned:
**     dl,db       double ecliptic longitude and latitude (radians)
**
**  1) The TT date date1+date2 is a Julian Date, apportioned in any
**     convenient way between the two arguments.  For example,
**     JD(TT)=2450123.7 could be expressed in any of these ways,
**     among others:
**
**            date1          date2
**
**         2450123.7           0.0       (JD method)
**         2451545.0       -1421.3       (J2000 method)
**         2400000.5       50123.2       (MJD method)
**         2450123.5           0.2       (date & time method)
**
**     The JD method is the most natural and convenient to use in
**     cases where the loss of several decimal digits of resolution
**     is acceptable.  The J2000 method is best matched to the way
**     the argument is handled internally and will deliver the
**     optimum resolution.  The MJD method and the date & time methods
**     are both good compromises between resolution and convenience.
**
**  2) No assumptions are made about whether the coordinates represent
**     starlight and embody astrometric effects such as parallax or
**     aberration.
**
**  3) The transformation is approximately that from mean J2000.0 right
**     ascension and declination to ecliptic longitude and latitude
**     (mean equinox and ecliptic of date), with only frame bias (always
**     less than 25 mas) to disturb this classical picture.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauEcm06     J2000.0 to ecliptic rotation matrix, IAU 2006
**     iauRxp       product of r-matrix and p-vector
**     iauC2s       unit vector to spherical coordinates
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dl = c_double()
    db = c_double()
    _sofa.iauEqec06(date1,date2,dr,dd,byref(dl),byref(db))
    return dl.value, db.value


#iauLteceq
try:
    _sofa.iauLteceq.argtypes = [
        c_double, # epj
        c_double, # dl
        c_double, # db
        POINTER(c_double), # dr
        POINTER(c_double), # dd
        ]
except AttributeError:
    pass
def lteceq(epj,dl,db):
    """
    Transformation from ecliptic coordinates (mean equinox and ecliptic  of date) to ICRS RA,Dec, using a long-term precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT).
    dl : float
        ecliptic longitude and latitude (radians).
    db : float
        ecliptic longitude and latitude (radians).

    Returns
    -------
   dr :  float
        ICRS right ascension  (radians).
    dd : float
        ICRS  declination (radians).
/*
**  - - - - - - - - - -
**   i a u L t e c e q
**  - - - - - - - - - -
**
**  Transformation from ecliptic coordinates (mean equinox and ecliptic
**  of date) to ICRS RA,Dec, using a long-term precession model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double     Julian epoch (TT)
**     dl,db   double     ecliptic longitude and latitude (radians)
**
**  Returned:
**     dr,dd   double     ICRS right ascension and declination (radians)
**
**  1) No assumptions are made about whether the coordinates represent
**     starlight and embody astrometric effects such as parallax or
**     aberration.
**
**  2) The transformation is approximately that from ecliptic longitude
**     and latitude (mean equinox and ecliptic of date) to mean J2000.0
**     right ascension and declination, with only frame bias (always
**     less than 25 mas) to disturb this classical picture.
**
**  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauLtecm     J2000.0 to ecliptic rotation matrix, long term
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauC2s       unit vector to spherical coordinates
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dr = c_double()
    dd = c_double()
    _sofa.iauLteceq(epj,dl,db,byref(dr),byref(dd))
    return dr.value, dd.value


#iauLtecm
try:
    _sofa.iauLtecm.argtypes = [
        c_double, # epj
        ndpointer(shape=(3,3),dtype = float),#rm[3][3]
        ]
except AttributeError:
    pass
def ltecm(epj):
    """
    ICRS equatorial to ecliptic rotation matrix, long-term.

    Parameters
    ----------
    epj : float
        DESCRIPTION.

    Returns
    -------
    rm : double[3][3]
        ICRS to ecliptic rotation matrix.
/*
**  - - - - - - - - -
**   i a u L t e c m
**  - - - - - - - - -
**
**  ICRS equatorial to ecliptic rotation matrix, long-term.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double         Julian epoch (TT)
**
**  Returned:
**     rm      double[3][3]   ICRS to ecliptic rotation matrix
**
**  Notes:
**
**  1) The matrix is in the sense
**
**        E_ep = rm x P_ICRS,
**
**     where P_ICRS is a vector with respect to ICRS right ascension
**     and declination axes and E_ep is the same vector with respect to
**     the (inertial) ecliptic and equinox of epoch epj.
**
**  2) P_ICRS is a free vector, merely a direction, typically of unit
**     magnitude, and not bound to any particular spatial origin, such
**     as the Earth, Sun or SSB.  No assumptions are made about whether
**     it represents starlight and embodies astrometric effects such as
**     parallax or aberration.  The transformation is approximately that
**     between mean J2000.0 right ascension and declination and ecliptic
**     longitude and latitude, with only frame bias (always less than
**     25 mas) to disturb this classical picture.
**
**  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  Called:
**     iauLtpequ    equator pole, long term
**     iauLtpecl    ecliptic pole, long term
**     iauPxp       vector product
**     iauPn        normalize vector
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2015 December 6
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rm = zeros(shape=(3,3),dtype=float)
    _sofa.iauLtecm(epj,rm)
    return rm


#iauLteqec
try:
    _sofa.iauLteqec.argtypes = [
        c_double, # epj
        c_double, # dr
        c_double, # dd
        POINTER(c_double), # dl
        POINTER(c_double), # db
        ]
except AttributeError:
    pass
def lteqec(epj,dr,dd):
    """
    Transformation from ICRS equatorial coordinates to ecliptic coordinates (mean equinox and ecliptic of date) using a long-term  precession model.

    Parameters
    ----------
    epj : float
        Julian epoch (TT).
    dr : float
        ICRS right ascension and declination (radians).
    dd : float
        ICRS right ascension and declination (radians).

    Returns
    -------
    dl : float
        ecliptic longitude and latitude (radians).
    db : float
        ecliptic longitude and latitude (radians).
/*
**  - - - - - - - - - -
**   i a u L t e q e c
**  - - - - - - - - - -
**
**  Transformation from ICRS equatorial coordinates to ecliptic
**  coordinates (mean equinox and ecliptic of date) using a long-term
**  precession model.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     epj     double     Julian epoch (TT)
**     dr,dd   double     ICRS right ascension and declination (radians)
**
**  Returned:
**     dl,db   double     ecliptic longitude and latitude (radians)
**
**  1) No assumptions are made about whether the coordinates represent
**     starlight and embody astrometric effects such as parallax or
**     aberration.
**
**  2) The transformation is approximately that from mean J2000.0 right
**     ascension and declination to ecliptic longitude and latitude
**     (mean equinox and ecliptic of date), with only frame bias (always
**     less than 25 mas) to disturb this classical picture.
**
**  3) The Vondrak et al. (2011, 2012) 400 millennia precession model
**     agrees with the IAU 2006 precession at J2000.0 and stays within
**     100 microarcseconds during the 20th and 21st centuries.  It is
**     accurate to a few arcseconds throughout the historical period,
**     worsening to a few tenths of a degree at the end of the
**     +/- 200,000 year time span.
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauLtecm     J2000.0 to ecliptic rotation matrix, long term
**     iauRxp       product of r-matrix and p-vector
**     iauC2s       unit vector to spherical coordinates
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**
**  References:
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2011, New precession
**    expressions, valid for long time intervals, Astron.Astrophys. 534,
**    A22
**
**    Vondrak, J., Capitaine, N. and Wallace, P., 2012, New precession
**    expressions, valid for long time intervals (Corrigendum),
**    Astron.Astrophys. 541, C1
**
**  This revision:  2016 February 9
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dl = c_double()
    db = c_double()
    _sofa.iauLteqec(epj,dr,dd,byref(dl),byref(db))
    return dl.value, db.value



#########################################################################
##   Galactic coordinates
#######################################################################
#iauG2icrs
try:
    _sofa.iauG2icrs.argtypes = [
        c_double, # dl
        c_double, # db
        POINTER(c_double), # dr
        POINTER(c_double), # dd
        ]
except AttributeError:
    pass
def g2icrs(dl,db):
    """
    Transformation from Galactic Coordinates to ICRS.

    Parameters
    ----------
    dl : float
        galactic longitude (radians).
    db : float
        galactic latitude (radians).

    Returns
    -------
    dr : float
        ICRS right ascension (radians).
    dd : float
        ICRS declination (radians).
/*
**  - - - - - - - - - -
**   i a u G 2 i c r s
**  - - - - - - - - - -
**
**  Transformation from Galactic Coordinates to ICRS.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     dl     double      galactic longitude (radians)
**     db     double      galactic latitude (radians)
**
**  Returned:
**     dr     double      ICRS right ascension (radians)
**     dd     double      ICRS declination (radians)
**
**  Notes:
**
**  1) The IAU 1958 system of Galactic coordinates was defined with
**     respect to the now obsolete reference system FK4 B1950.0.  When
**     interpreting the system in a modern context, several factors have
**     to be taken into account:
**
**     . The inclusion in FK4 positions of the E-terms of aberration.
**
**     . The distortion of the FK4 proper motion system by differential
**       Galactic rotation.
**
**     . The use of the B1950.0 equinox rather than the now-standard
**       J2000.0.
**
**     . The frame bias between ICRS and the J2000.0 mean place system.
**
**     The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
**     matrix that transforms directly between ICRS and Galactic
**     coordinates with the above factors taken into account.  The
**     matrix is derived from three angles, namely the ICRS coordinates
**     of the Galactic pole and the longitude of the ascending node of
**     the galactic equator on the ICRS equator.  They are given in
**     degrees to five decimal places and for canonical purposes are
**     regarded as exact.  In the Hipparcos Catalogue the matrix
**     elements are given to 10 decimal places (about 20 microarcsec).
**     In the present SOFA function the matrix elements have been
**     recomputed from the canonical three angles and are given to 30
**     decimal places.
**
**  2) The inverse transformation is performed by the function iauIcrs2g.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**     iauS2c       spherical coordinates to unit vector
**     iauTrxp      product of transpose of r-matrix and p-vector
**     iauC2s       p-vector to spherical
**
**  Reference:
**     Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
**     catalogues.  Astrometric and photometric star catalogues
**     derived from the ESA Hipparcos Space Astrometry Mission.  ESA
**     Publications Division, Noordwijk, Netherlands.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dr = c_double()
    dd = c_double()
    _sofa.iauG2icrs(dl,db,byref(dr),byref(dd))
    return dr.value, dd.value


#iauIcrs2g
try:
    _sofa.iauIcrs2g.argtypes = [
        c_double, # dr
        c_double, # dd
        POINTER(c_double), # dl
        POINTER(c_double), # db
        ]
except AttributeError:
    pass
def icrs2g(dr,dd):
    """
    Transformation from ICRS to Galactic Coordinates.

    Parameters
    ----------
    dr : float
        ICRS right ascension (radians).
    dd : float
        ICRS declination (radians).

    Returns
    -------
    dl : float
        galactic longitude (radians).
    db : float
        galactic latitude (radians).
/*
**  - - - - - - - - - -
**   i a u I c r s 2 g
**  - - - - - - - - - -
**
**  Transformation from ICRS to Galactic Coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     dr     double      ICRS right ascension (radians)
**     dd     double      ICRS declination (radians)
**
**  Returned:
**     dl     double      galactic longitude (radians)
**     db     double      galactic latitude (radians)
**
**  Notes:
**
**  1) The IAU 1958 system of Galactic coordinates was defined with
**     respect to the now obsolete reference system FK4 B1950.0.  When
**     interpreting the system in a modern context, several factors have
**     to be taken into account:
**
**     . The inclusion in FK4 positions of the E-terms of aberration.
**
**     . The distortion of the FK4 proper motion system by differential
**       Galactic rotation.
**
**     . The use of the B1950.0 equinox rather than the now-standard
**       J2000.0.
**
**     . The frame bias between ICRS and the J2000.0 mean place system.
**
**     The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
**     matrix that transforms directly between ICRS and Galactic
**     coordinates with the above factors taken into account.  The
**     matrix is derived from three angles, namely the ICRS coordinates
**     of the Galactic pole and the longitude of the ascending node of
**     the galactic equator on the ICRS equator.  They are given in
**     degrees to five decimal places and for canonical purposes are
**     regarded as exact.  In the Hipparcos Catalogue the matrix
**     elements are given to 10 decimal places (about 20 microarcsec).
**     In the present SOFA function the matrix elements have been
**     recomputed from the canonical three angles and are given to 30
**     decimal places.
**
**  2) The inverse transformation is performed by the function iauG2icrs.
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**     iauAnpm      normalize angle into range +/- pi
**     iauS2c       spherical coordinates to unit vector
**     iauRxp       product of r-matrix and p-vector
**     iauC2s       p-vector to spherical
**
**  Reference:
**     Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
**     catalogues.  Astrometric and photometric star catalogues
**     derived from the ESA Hipparcos Space Astrometry Mission.  ESA
**     Publications Division, Noordwijk, Netherlands.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    dl = c_double()
    db = c_double()
    _sofa.iauIcrs2g(dr,dd,byref(dl),byref(db))
    return dl.value, db.value



#########################################################################
###   Geodetic/geocentric
####################################################################
  
#iauEform
eform_msg = {0: 'OK', # Unused
                -1:'illegal identifier (Note 3)',}
  
try:
    _sofa.iauEform.argtypes = [
    c_int, # n
    POINTER(c_double), # a
    POINTER(c_double), # f
    ]
except AttributeError:
    pass
def eform(n):
    """
    Earth reference ellipsoids.

    Parameters
    ----------
    n : int
        ellipsoid identifier (Note 1).

    Returns
    -------
    a : float
        DEequatorial radius (meters, Note 2)SCRIPTION.
    f : float
        flattening (Note 2).
/*
**  - - - - - - - - -
**   i a u E f o r m
**  - - - - - - - - -
**
**  Earth reference ellipsoids.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical.
**
**  Given:
**     n    int         ellipsoid identifier (Note 1)
**
**  Returned:
**     a    double      equatorial radius (meters, Note 2)
**     f    double      flattening (Note 2)
**
**  Returned (function value):
**          int         status:  0 = OK
**                              -1 = illegal identifier (Note 3)
**
**  Notes:
**
**  1) The identifier n is a number that specifies the choice of
**     reference ellipsoid.  The following are supported:
**
**        n    ellipsoid
**
**        1     WGS84
**        2     GRS80
**        3     WGS72
**
**     The n value has no significance outside the SOFA software.  For
**     convenience, symbols WGS84 etc. are defined in sofam.h.
**
**  2) The ellipsoid parameters are returned in the form of equatorial
**     radius in meters (a) and flattening (f).  The latter is a number
**     around 0.00335, i.e. around 1/298.
**
**  3) For the case where an unsupported n value is supplied, zero a and
**     f are returned, as well as error status.
**
**  References:
**
**     Department of Defense World Geodetic System 1984, National
**     Imagery and Mapping Agency Technical Report 8350.2, Third
**     Edition, p3-2.
**
**     Moritz, H., Bull. Geodesique 66-2, 187 (1992).
**
**     The Department of Defense World Geodetic System 1972, World
**     Geodetic System Committee, May 1974.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     p220.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = c_double()
    f = c_double()
    s = _sofa.iauEform(n,byref(a),byref(f))
    if s != 0:
        warnings.warn(eform_msg[s], UserWarning, 2)
    return a.value, f.value


#iauGc2gd
gc2gd_msg = {0: 'OK', # Unused
                -1:'illegal identifier (Note 3)',
                -2 : 'internal error (Note 3)'}
try:
    _sofa.iauGc2gd.argtypes = [
    c_int, # n
    ndpointer(shape = 3,dtype = float),#xyz[3]
    POINTER(c_double), # elong
    POINTER(c_double), # phi
    POINTER(c_double),# height
    ]
except AttributeError:
    pass
def gc2gd(n,xyz):
    """
    Transform geocentric coordinates to geodetic using the specified reference ellipsoid

    Parameters
    ----------
    n : int
        ellipsoid identifier (Note 1).
    xyz : double[3]
        geocentric vector (Note 2).

    Returns
    -------
    elong : float
        longitude (radians, east +ve, Note 3).
    phi : float
        latitude (geodetic, radians, Note 3).
    height : float
        height above ellipsoid (geodetic, Notes 2,3).
/*
**  - - - - - - - - -
**   i a u G c 2 g d
**  - - - - - - - - -
**
**  Transform geocentric coordinates to geodetic using the specified
**  reference ellipsoid.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical transformation.
**
**  Given:
**     n       int        ellipsoid identifier (Note 1)
**     xyz     double[3]  geocentric vector (Note 2)
**
**  Returned:
**     elong   double     longitude (radians, east +ve, Note 3)
**     phi     double     latitude (geodetic, radians, Note 3)
**     height  double     height above ellipsoid (geodetic, Notes 2,3)
**
**  Returned (function value):
**            int         status:  0 = OK
**                                -1 = illegal identifier (Note 3)
**                                -2 = internal error (Note 3)
**
**  Notes:
**
**  1) The identifier n is a number that specifies the choice of
**     reference ellipsoid.  The following are supported:
**
**        n    ellipsoid
**
**        1     WGS84
**        2     GRS80
**        3     WGS72
**
**     The n value has no significance outside the SOFA software.  For
**     convenience, symbols WGS84 etc. are defined in sofam.h.
**
**  2) The geocentric vector (xyz, given) and height (height, returned)
**     are in meters.
**
**  3) An error status -1 means that the identifier n is illegal.  An
**     error status -2 is theoretically impossible.  In all error cases,
**     all three results are set to -1e9.
**
**  4) The inverse transformation is performed in the function iauGd2gc.
**
**  Called:
**     iauEform     Earth reference ellipsoids
**     iauGc2gde    geocentric to geodetic transformation, general
**
**  This revision:  2013 September 1
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    elong = c_double()
    phi = c_double()
    height = c_double()
    s = _sofa.iauGc2gd(n,xyz,byref(elong),byref(phi),byref(height))
    if s != 0:
        warnings.warn(gc2gd_msg[s], UserWarning, 2)
    return elong.value, phi.value, height.value


#iauGc2gde
gc2gde_msg = {0: 'OK', # Unused
                -1:'illegal f',
                -2 : 'illegal a'}
try:
    _sofa.iauGc2gde.argtypes = [
    c_double, # a
    c_double, # f
    ndpointer(shape = 3,dtype = float), # xyz[3]
    POINTER(c_double), # elong
    POINTER(c_double), # phi
    POINTER(c_double),# height
    ]
except AttributeError:
    pass
def gc2gde(a,f,xyz):
    """
    Transform geocentric coordinates to geodetic for a reference  ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius (Notes 2,4).
    f : float
        flattening (Note 3).
    xyz : double[3]
        geocentric vector (Note 4).

    Returns
    -------
    elong : float
        longitude (radians, east +ve).
    phi : float
        latitude (geodetic, radians).
    height : float
        height above ellipsoid (geodetic, Note 4).
/*
**  - - - - - - - - - -
**   i a u G c 2 g d e
**  - - - - - - - - - -
**
**  Transform geocentric coordinates to geodetic for a reference
**  ellipsoid of specified form.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     a       double     equatorial radius (Notes 2,4)
**     f       double     flattening (Note 3)
**     xyz     double[3]  geocentric vector (Note 4)
**
**  Returned:
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians)
**     height  double     height above ellipsoid (geodetic, Note 4)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal f
**                                -2 = illegal a
**
**  Notes:
**
**  1) This function is based on the GCONV2H Fortran subroutine by
**     Toshio Fukushima (see reference).
**
**  2) The equatorial radius, a, can be in any units, but meters is
**     the conventional choice.
**
**  3) The flattening, f, is (for the Earth) a value around 0.00335,
**     i.e. around 1/298.
**
**  4) The equatorial radius, a, and the geocentric vector, xyz,
**     must be given in the same units, and determine the units of
**     the returned height, height.
**
**  5) If an error occurs (status < 0), elong, phi and height are
**     unchanged.
**
**  6) The inverse transformation is performed in the function
**     iauGd2gce.
**
**  7) The transformation for a standard ellipsoid (such as WGS84) can
**     more conveniently be performed by calling iauGc2gd, which uses a
**     numerical code to identify the required A and F values.
**
**  Reference:
**
**     Fukushima, T., "Transformation from Cartesian to geodetic
**     coordinates accelerated by Halley's method", J.Geodesy (2006)
**     79: 689-693
**
**  This revision:  2014 November 7
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    elong = c_double()
    phi = c_double()
    height = c_double()
    s = _sofa.iauGc2gde(a,f,xyz,byref(elong),byref(phi),byref(height))
    if s != 0:
        warnings.warn(gc2gde_msg[s], UserWarning, 2)
    return elong.value, phi.value, height.value
 

#iauGd2gc
gd2gc_msg = {0: 'OK', # Unused
                -1:'illegal identifier (Note 3)',
                -2 : 'illegal case (Note 3)'}
try:
    _sofa.iauGd2gc.argtypes = [
    c_int, # n
    c_double, # elong
    c_double,# phi
    c_double,# height
    ndpointer(shape = 3,dtype = float),#xyz[3]
    ]
except AttributeError:
    pass
def gd2gc(n,elong,phi,height):
    """
    Transform geodetic coordinates to geocentric using the specified reference ellipsoid.

    Parameters
    ----------
    n : int
        ellipsoid identifier (Note 1).
    elong : float
        longitude (radians, east +ve).
    phi : float
        latitude (geodetic, radians, Note 3).
    height : float
        height above ellipsoid (geodetic, Notes 2,3).

    Returns
    -------
    xyz : double[3]
        geocentric vector (Note 2).
/*
**  - - - - - - - - -
**   i a u G d 2 g c
**  - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric using the specified
**  reference ellipsoid.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  canonical transformation.
**
**  Given:
**     n       int        ellipsoid identifier (Note 1)
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians, Note 3)
**     height  double     height above ellipsoid (geodetic, Notes 2,3)
**
**  Returned:
**     xyz     double[3]  geocentric vector (Note 2)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal identifier (Note 3)
**                                -2 = illegal case (Note 3)
**
**  Notes:
**
**  1) The identifier n is a number that specifies the choice of
**     reference ellipsoid.  The following are supported:
**
**        n    ellipsoid
**
**        1     WGS84
**        2     GRS80
**        3     WGS72
**
**     The n value has no significance outside the SOFA software.  For
**     convenience, symbols WGS84 etc. are defined in sofam.h.
**
**  2) The height (height, given) and the geocentric vector (xyz,
**     returned) are in meters.
**
**  3) No validation is performed on the arguments elong, phi and
**     height.  An error status -1 means that the identifier n is
**     illegal.  An error status -2 protects against cases that would
**     lead to arithmetic exceptions.  In all error cases, xyz is set
**     to zeros.
**
**  4) The inverse transformation is performed in the function iauGc2gd.
**
**  Called:
**     iauEform     Earth reference ellipsoids
**     iauGd2gce    geodetic to geocentric transformation, general
**     iauZp        zero p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    xyz = zeros(shape = 3,dtype = float)
    s = _sofa.iauGd2gc(n,elong,phi,height,xyz)
    if s != 0:
        warnings.warn(gd2gc_msg[s], UserWarning, 2)
    return xyz


#iauGd2gce
gd2gce_msg = {0: 'OK', # Unused
                -1:'illegal case (Note 4)'}
try:
    _sofa.iauGd2gce.argtypes = [
    c_double, # a
    c_double, # f
    c_double, # elong
    c_double, # phi
    c_double,# height
    ndpointer(shape = 3,dtype = float),#xyz[3]
    ]
except AttributeError:
    pass
def gd2gce(a,f,elong,phi,height):
    """
    Transform geodetic coordinates to geocentric for a reference ellipsoid of specified form.

    Parameters
    ----------
    a : float
        equatorial radius (Notes 1,4).
    f : float
        flattening (Notes 2,4).
    elong : float
        longitude (radians, east +ve).
    phi : float
        latitude (geodetic, radians, Note 4).
    height : float
        height above ellipsoid (geodetic, Notes 3,4).

    Returns
    -------
    xyz : double[3]
        geocentric vector (Note 3).
/*
**  - - - - - - - - - -
**   i a u G d 2 g c e
**  - - - - - - - - - -
**
**  Transform geodetic coordinates to geocentric for a reference
**  ellipsoid of specified form.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     a       double     equatorial radius (Notes 1,4)
**     f       double     flattening (Notes 2,4)
**     elong   double     longitude (radians, east +ve)
**     phi     double     latitude (geodetic, radians, Note 4)
**     height  double     height above ellipsoid (geodetic, Notes 3,4)
**
**  Returned:
**     xyz     double[3]  geocentric vector (Note 3)
**
**  Returned (function value):
**             int        status:  0 = OK
**                                -1 = illegal case (Note 4)
**  Notes:
**
**  1) The equatorial radius, a, can be in any units, but meters is
**     the conventional choice.
**
**  2) The flattening, f, is (for the Earth) a value around 0.00335,
**     i.e. around 1/298.
**
**  3) The equatorial radius, a, and the height, height, must be
**     given in the same units, and determine the units of the
**     returned geocentric vector, xyz.
**
**  4) No validation is performed on individual arguments.  The error
**     status -1 protects against (unrealistic) cases that would lead
**     to arithmetic exceptions.  If an error occurs, xyz is unchanged.
**
**  5) The inverse transformation is performed in the function
**     iauGc2gde.
**
**  6) The transformation for a standard ellipsoid (such as WGS84) can
**     more conveniently be performed by calling iauGd2gc,  which uses a
**     numerical code to identify the required a and f values.
**
**  References:
**
**     Green, R.M., Spherical Astronomy, Cambridge University Press,
**     (1985) Section 4.5, p96.
**
**     Explanatory Supplement to the Astronomical Almanac,
**     P. Kenneth Seidelmann (ed), University Science Books (1992),
**     Section 4.22, p202.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    xyz = zeros(shape = 3,dtype = float)
    s = _sofa.iauGd2gce(a,f,elong,phi,height,xyz)
    if s != 0:
        warnings.warn(gd2gce_msg[s], UserWarning, 2)
    return xyz



#########################################################################
###   Gnomonic projection
####################################################################
    
#iauTpors
tpors_msg = {0: 'OK', # Unused
                1:'only the first solution is useful (Note 6)',
                2 : 'both solutions are useful (Note 6)'}
try:
    _sofa.iauTpors.argtypes = [
    c_double, # xi
    c_double, # eta
    c_double, # a
    c_double, # b
    POINTER(c_double), # a01
    POINTER(c_double), # b01
    POINTER(c_double),# a02
    POINTER(c_double),# b02
    ]
except AttributeError:
    pass
def tpors(xi,eta,a,b):
    """
    In the tangent plane projection, given the rectangular coordinates of a star and its spherical coordinates, determine the spherical  coordinates of the tangent point.

    Parameters
    ----------
    xi : double
        rectangular coordinates of star image (Note 2).
    eta : double
        rectangular coordinates of star image (Note 2).
    a : double
        star's spherical coordinates (Note 3).
    b : double
        star's spherical coordinates (Note 3).

    Returns
    -------
    a01 : double
        tangent point's spherical coordinates, Soln. 1.
    b01 :double
        tangent point's spherical coordinates, Soln. 1.
    double
        tangent point's spherical coordinates, Soln. 2.
    double
        tangent point's spherical coordinates, Soln. 2.
/*
**  - - - - - - - - -
**   i a u T p o r s
**  - - - - - - - - -
**
**  In the tangent plane projection, given the rectangular coordinates
**  of a star and its spherical coordinates, determine the spherical
**  coordinates of the tangent point.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xi,eta     double  rectangular coordinates of star image (Note 2)
**     a,b        double  star's spherical coordinates (Note 3)
**
**  Returned:
**     *a01,*b01  double  tangent point's spherical coordinates, Soln. 1
**     *a02,*b02  double  tangent point's spherical coordinates, Soln. 2
**
**  Returned (function value):
**                int     number of solutions:
**                        0 = no solutions returned (Note 5)
**                        1 = only the first solution is useful (Note 6)
**                        2 = both solutions are useful (Note 6)
**
**  Notes:
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the spherical coordinates are observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the spherical coordinates are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) All angular arguments are in radians.
**
**  4) The angles a01 and a02 are returned in the range 0-2pi.  The
**     angles b01 and b02 are returned in the range +/-pi, but in the
**     usual, non-pole-crossing, case, the range is +/-pi/2.
**
**  5) Cases where there is no solution can arise only near the poles.
**     For example, it is clearly impossible for a star at the pole
**     itself to have a non-zero xi value, and hence it is meaningless
**     to ask where the tangent point would have to be to bring about
**     this combination of xi and dec.
**
**  6) Also near the poles, cases can arise where there are two useful
**     solutions.  The return value indicates whether the second of the
**     two solutions returned is useful;  1 indicates only one useful
**     solution, the usual case.
**
**  7) The basis of the algorithm is to solve the spherical triangle PSC,
**     where P is the north celestial pole, S is the star and C is the
**     tangent point.  The spherical coordinates of the tangent point are
**     [a0,b0];  writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), side c
**     is then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
**     found) is (pi/2-b0).  Angle C is given by sin(C) = xi/rho and
**     cos(C) = eta/rho.  Angle P (to be found) is the longitude
**     difference between star and tangent point (a-a0).
**
**  8) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes      iauTpxev         xi,eta
**         iauTpsts      iauTpstv          star
**       > iauTpors <    iauTporv         origin
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a01 = c_double()
    b01 = c_double()
    a02 = c_double()
    b02 = c_double()
    s = _sofa.iauTpors(xi,eta,a,b,byref(a01),byref(b01),byref(a02),byref(b02))
    if s != 0:
        warnings.warn(tpors_msg[s], UserWarning, 2)
    return a01.value, b01.value, a02.value, b02.value


#iauTporv
tporv_msg = {0: 'OK', # Unused
                1:'only the first solution is useful (Note 5)',
                2 : 'both solutions are useful (Note 5)'}
try:
    _sofa.iauTporv.argtypes = [
    c_double, # xi
    c_double, # eta
    ndpointer(shape = 3, dtype = float), #v[3]
    ndpointer(shape = 3, dtype = float), #v01[3]
    ndpointer(shape = 3, dtype = float), #v02[3]
    ]
except AttributeError:
    pass
def tporv(xi,eta,v):
    """
    In the tangent plane projection, given the rectangular coordinates of a star and its direction cosines, determine the direction cosines of the tangent point.

    Parameters
    ----------
    xi : double
        rectangular coordinates of star image (Note 2).
    eta : double
        rectangular coordinates of star image (Note 2).
    v : double[3]
        star's direction cosines (Note 3).

    Returns
    -------
    v01 : double[3]
        tangent point's direction cosines, Solution 1.
    v02 : double[3]
        tangent point's direction cosines, Solution 2.
/*
**  - - - - - - - - -
**   i a u T p o r v
**  - - - - - - - - -
**
**  In the tangent plane projection, given the rectangular coordinates
**  of a star and its direction cosines, determine the direction
**  cosines of the tangent point.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xi,eta   double    rectangular coordinates of star image (Note 2)
**     v        double[3] star's direction cosines (Note 3)
**
**  Returned:
**     v01      double[3] tangent point's direction cosines, Solution 1
**     v02      double[3] tangent point's direction cosines, Solution 2
**
**  Returned (function value):
**                int     number of solutions:
**                        0 = no solutions returned (Note 4)
**                        1 = only the first solution is useful (Note 5)
**                        2 = both solutions are useful (Note 5)
**
**  Notes:
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the direction cosines represent observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the direction cosines are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) The vector v must be of unit length or the result will be wrong.
**
**  4) Cases where there is no solution can arise only near the poles.
**     For example, it is clearly impossible for a star at the pole
**     itself to have a non-zero xi value, and hence it is meaningless
**     to ask where the tangent point would have to be.
**
**  5) Also near the poles, cases can arise where there are two useful
**     solutions.  The return value indicates whether the second of the
**     two solutions returned is useful;  1 indicates only one useful
**     solution, the usual case.
**
**  6) The basis of the algorithm is to solve the spherical triangle
**     PSC, where P is the north celestial pole, S is the star and C is
**     the tangent point.  Calling the celestial spherical coordinates
**     of the star and tangent point (a,b) and (a0,b0) respectively, and
**     writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), and
**     transforming the vector v into (a,b) in the normal way, side c is
**     then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
**     found) is (pi/2-b0), while angle C is given by sin(C) = xi/rho
**     and cos(C) = eta/rho;  angle P (to be found) is (a-a0).  After
**     solving the spherical triangle, the result (a0,b0) can be
**     expressed in vector form as v0.
**
**  7) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes      iauTpxev         xi,eta
**         iauTpsts      iauTpstv          star
**         iauTpors    > iauTporv <       origin
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    v01 = zeros(shape = 3, dtype = float)
    v02 = zeros(shape = 3, dtype = float)
    s = _sofa.iauTporv(xi,eta,v,v01,v02)
    if s != 0:
        warnings.warn(tporv_msg[s], UserWarning, 2)
    return v01,v02


#iauTpsts
try:
    _sofa.iauTpsts.argtypes = [
    c_double, # xi
    c_double, # eta
    c_double, # a0
    c_double, # b0
    POINTER(c_double),# a
    POINTER(c_double),# b
    ]
except AttributeError:
    pass
def tpsts(xi,eta,a0,b0):
    """
    In the tangent plane projection, given the star's rectangular  coordinates and the spherical coordinates of the tangent point, solve for the spherical coordinates of the star.

    Parameters
    ----------
    xi : double
        rectangular coordinates of star image (Note 2).
    eta : double
        rectangular coordinates of star image (Note 2).
    a0 : double
        tangent point's spherical coordinates.
    b0 : double
        tangent point's spherical coordinates.

    Returns
    -------
    a : double
        star's spherical coordinates.
    b : double
        star's spherical coordinates.
/*
**  - - - - - - - - -
**   i a u T p s t s
**  - - - - - - - - -
**
**  In the tangent plane projection, given the star's rectangular
**  coordinates and the spherical coordinates of the tangent point,
**  solve for the spherical coordinates of the star.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xi,eta    double  rectangular coordinates of star image (Note 2)
**     a0,b0     double  tangent point's spherical coordinates
**
**  Returned:
**     *a,*b     double  star's spherical coordinates
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the spherical coordinates are observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the spherical coordinates are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) All angular arguments are in radians.
**
**  4) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes      iauTpxev         xi,eta
**       > iauTpsts <    iauTpstv          star
**         iauTpors      iauTporv         origin
**
**  Called:
**     iauAnp       normalize angle into range 0 to 2pi
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    a = c_double()
    b = c_double()
    _sofa.iauTpsts(xi,eta,a0,b0,byref(a),byref(b))
    return  a.value, b.value


#iauTpstv
try:
    _sofa.iauTpstv.argtypes = [
    c_double, # xi
    c_double, # eta
    ndpointer(shape=3,dtype =float),#v0[3]
    ndpointer(shape=3,dtype =float),#v[3]
    ]
except AttributeError:
    pass
def tpstv(xi,eta,v0):
    """
    In the tangent plane projection, given the star's rectangular  coordinates and the direction cosines of the tangent point, solve  for the direction cosines of the star.

    Parameters
    ----------
    xi : double
        rectangular coordinates of star image (Note 2).
    eta : double
        rectangular coordinates of star image (Note 2).
    v0 : double[3]
        tangent point's direction cosines.

    Returns
    -------
    v : double[3]
        star's direction cosines.
/*
**  - - - - - - - - -
**   i a u T p s t v
**  - - - - - - - - -
**
**  In the tangent plane projection, given the star's rectangular
**  coordinates and the direction cosines of the tangent point, solve
**  for the direction cosines of the star.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     xi,eta  double     rectangular coordinates of star image (Note 2)
**     v0      double[3]  tangent point's direction cosines
**
**  Returned:
**     v       double[3]  star's direction cosines
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the direction cosines represent observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the direction cosines are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) The method used is to complete the star vector in the (xi,eta)
**     based triad and normalize it, then rotate the triad to put the
**     tangent point at the pole with the x-axis aligned to zero
**     longitude.  Writing (a0,b0) for the celestial spherical
**     coordinates of the tangent point, the sequence of rotations is
**     (b-pi/2) around the x-axis followed by (-a-pi/2) around the
**     z-axis.
**
**  4) If vector v0 is not of unit length, the returned vector v will
**     be wrong.
**
**  5) If vector v0 points at a pole, the returned vector v will be
**     based on the arbitrary assumption that the longitude coordinate
**     of the tangent point is zero.
**
**  6) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes      iauTpxev         xi,eta
**         iauTpsts    > iauTpstv <        star
**         iauTpors      iauTporv         origin
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    v = zeros(shape=3,dtype = float)
    _sofa.iauTpstv(xi,eta,v0,v)
    return  v


#iauTpxes
tpxes_msg = {0 : 'OK', # Unused
                1 :'star too far from axis',
                2 : 'antistar on tangent plane',
                3 : 'antistar too far from axis'}
try:
    _sofa.iauTpxes.argtypes = [
    c_double, # a
    c_double, # b
    c_double, # a0
    c_double, # b0
    POINTER(c_double), # xi
    POINTER(c_double), # eta
    ]
except AttributeError:
    pass
def tpxes(a,b,a0,b0):
    """
    In the tangent plane projection, given celestial spherical  coordinates for a star and the tangent point, solve for the star's  rectangular coordinates in the tangent plane.

    Parameters
    ----------
    a : double
        star's spherical coordinates.
    b : double
        star's spherical coordinates.
    a0 : double
        tangent point's spherical coordinates.
    b0 : double
        tangent point's spherical coordinates.

    Returns
    -------
    double
        rectangular coordinates of star image (Note 2).
    double
        rectangular coordinates of star image (Note 2).
/*
**  - - - - - - - - -
**   i a u T p x e s
**  - - - - - - - - -
**
**  In the tangent plane projection, given celestial spherical
**  coordinates for a star and the tangent point, solve for the star's
**  rectangular coordinates in the tangent plane.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     a,b       double  star's spherical coordinates
**     a0,b0     double  tangent point's spherical coordinates
**
**  Returned:
**     *xi,*eta  double  rectangular coordinates of star image (Note 2)
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = star too far from axis
**                                2 = antistar on tangent plane
**                                3 = antistar too far from axis
**
**  Notes:
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the spherical coordinates are observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  For right-handed spherical coordinates,
**     (xi,eta) are also right-handed.  The units of (xi,eta) are,
**     effectively, radians at the tangent point.
**
**  3) All angular arguments are in radians.
**
**  4) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**       > iauTpxes <    iauTpxev         xi,eta
**         iauTpsts      iauTpstv          star
**         iauTpors      iauTporv         origin
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    xi = c_double()
    eta = c_double()
    s = _sofa.iauTpxes(a,b,a0,b0,byref(xi),byref(eta))
    if s != 0:
        warnings.warn(tpxes_msg[s], UserWarning, 2)
    return xi.value, eta.value


#iauTpxev
tpxev_msg = {0 : 'OK', # Unused
                1 :'star too far from axis',
                2 : 'antistar on tangent plane',
                3 : 'antistar too far from axis'}
try:
    _sofa.iauTpxev.argtypes = [
    ndpointer(shape=3,dtype=float), # v[3] 
    ndpointer(shape=3,dtype=float), # v0[3]
    POINTER(c_double), # xi
    POINTER(c_double), # eta
    ]
except AttributeError:
    pass
def tpxev(v,v0):
    """
    In the tangent plane projection, given celestial direction cosines  for a star and the tangent point, solve for the star's rectangular  coordinates in the tangent plane.

    Parameters
    ----------
    v : double[3]
        direction cosines of star (Note 4).
    v0 : double[3]
        direction cosines of tangent point (Note 4).

    Returns
    -------
    v : double
        tangent plane coordinates of star.
    v0 : double
        tangent plane coordinates of star.
/*
**  - - - - - - - - -
**   i a u T p x e v
**  - - - - - - - - -
**
**  In the tangent plane projection, given celestial direction cosines
**  for a star and the tangent point, solve for the star's rectangular
**  coordinates in the tangent plane.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     v         double[3]  direction cosines of star (Note 4)
**     v0        double[3]  direction cosines of tangent point (Note 4)
**
**  Returned:
**     *xi,*eta  double     tangent plane coordinates of star
**
**  Returned (function value):
**               int        status: 0 = OK
**                                  1 = star too far from axis
**                                  2 = antistar on tangent plane
**                                  3 = antistar too far from axis
**
**  Notes:
**
**  1) The tangent plane projection is also called the "gnomonic
**     projection" and the "central projection".
**
**  2) The eta axis points due north in the adopted coordinate system.
**     If the direction cosines represent observed (RA,Dec), the tangent
**     plane coordinates (xi,eta) are conventionally called the
**     "standard coordinates".  If the direction cosines are with
**     respect to a right-handed triad, (xi,eta) are also right-handed.
**     The units of (xi,eta) are, effectively, radians at the tangent
**     point.
**
**  3) The method used is to extend the star vector to the tangent
**     plane and then rotate the triad so that (x,y) becomes (xi,eta).
**     Writing (a,b) for the celestial spherical coordinates of the
**     star, the sequence of rotations is (a+pi/2) around the z-axis
**     followed by (pi/2-b) around the x-axis.
**
**  4) If vector v0 is not of unit length, or if vector v is of zero
**     length, the results will be wrong.
**
**  5) If v0 points at a pole, the returned (xi,eta) will be based on
**     the arbitrary assumption that the longitude coordinate of the
**     tangent point is zero.
**
**  6) This function is a member of the following set:
**
**         spherical      vector         solve for
**
**         iauTpxes    > iauTpxev <       xi,eta
**         iauTpsts      iauTpstv          star
**         iauTpors      iauTporv         origin
**
**  References:
**
**     Calabretta M.R. & Greisen, E.W., 2002, "Representations of
**     celestial coordinates in FITS", Astron.Astrophys. 395, 1077
**
**     Green, R.M., "Spherical Astronomy", Cambridge University Press,
**     1987, Chapter 13.
**
**  This revision:   2018 January 2
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    xi = c_double()
    eta = c_double()
    s = _sofa.iauTpxev(v,v0,byref(xi),byref(eta))
    if s != 0:
        warnings.warn(tpxev_msg[s], UserWarning, 2)
    return xi.value, eta.value


#################################################################
###########     Horizon/equatoria
######################################################

#iauAe2hd
try:
    _sofa.iauAe2hd.argtypes = [
    c_double, # az
    c_double, # el
    c_double, # phi
    POINTER(c_double), # ha
    POINTER(c_double), # dec
    ]
except AttributeError:
    pass
def ae2hd(az,el,phi):
    """
    Horizon to equatorial coordinates:  transform azimuth and altitude  to hour angle and declination.

    Parameters
    ----------
    az : double
        azimuth.
    el : double
        altitude (informally, elevation).
    phi : double
        site latitude.

    Returns
    -------
    ha : double
        hour angle (local).
    dec : double
        declination.
/*
**  - - - - - - - - -
**   i a u A e 2 h d
**  - - - - - - - - -
**
**  Horizon to equatorial coordinates:  transform azimuth and altitude
**  to hour angle and declination.
**
**  Given:
**     az       double       azimuth
**     el       double       altitude (informally, elevation)
**     phi      double       site latitude
**
**  Returned:
**     ha       double       hour angle (local)
**     dec      double       declination
**
**  Notes:
**
**  1)  All the arguments are angles in radians.
**
**  2)  The sign convention for azimuth is north zero, east +pi/2.
**
**  3)  HA is returned in the range +/-pi.  Declination is returned in
**      the range +/-pi/2.
**
**  4)  The latitude phi is pi/2 minus the angle between the Earth's
**      rotation axis and the adopted zenith.  In many applications it
**      will be sufficient to use the published geodetic latitude of the
**      site.  In very precise (sub-arcsecond) applications, phi can be
**      corrected for polar motion.
**
**  5)  The azimuth az must be with respect to the rotational north pole,
**      as opposed to the ITRS pole, and an azimuth with respect to north
**      on a map of the Earth's surface will need to be adjusted for
**      polar motion if sub-arcsecond accuracy is required.
**
**  6)  Should the user wish to work with respect to the astronomical
**      zenith rather than the geodetic zenith, phi will need to be
**      adjusted for deflection of the vertical (often tens of
**      arcseconds), and the zero point of ha will also be affected.
**
**  7)  The transformation is the same as Ve = Ry(phi-pi/2)*Rz(pi)*Vh,
**      where Ve and Vh are lefthanded unit vectors in the (ha,dec) and
**      (az,el) systems respectively and Rz and Ry are rotations about
**      first the z-axis and then the y-axis.  (n.b. Rz(pi) simply
**      reverses the signs of the x and y components.)  For efficiency,
**      the algorithm is written out rather than calling other utility
**      functions.  For applications that require even greater
**      efficiency, additional savings are possible if constant terms
**      such as functions of latitude are computed once and for all.
**
**  8)  Again for efficiency, no range checking of arguments is carried
**      out.
**
**  Last revision:   2017 September 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ha = c_double()
    dec = c_double()
    _sofa.iauAe2hd(az,el,phi,byref(ha),byref(dec))
    return ha.value, dec.value

#iauHd2ae
try:
    _sofa.iauHd2ae.argtypes = [
    c_double, # ha
    c_double, # dec
    c_double, # phi
    POINTER(c_double), # az
    POINTER(c_double), # el
    ]
except AttributeError:
    pass
def hd2ae(ha, dec, phi):
    """
    Equatorial to horizon coordinates:  transform hour angle and  declination to azimuth and altitude.

    Parameters
    ----------
    ha : double
        hour angle (local).
    dec : double
        declination.
    phi : double
        site latitude.

    Returns
    -------
    az : double
        azimuth.
    el : double
        altitude(informally,elevation).
/*
**  - - - - - - - - -
**   i a u H d 2 a e
**  - - - - - - - - -
**
**  Equatorial to horizon coordinates:  transform hour angle and
**  declination to azimuth and altitude.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     ha       double       hour angle (local)
**     dec      double       declination
**     phi      double       site latitude
**
**  Returned:
**     *az      double       azimuth
**     *el      double       altitude (informally, elevation)
**
**  Notes:
**
**  1)  All the arguments are angles in radians.
**
**  2)  Azimuth is returned in the range 0-2pi;  north is zero, and east
**      is +pi/2.  Altitude is returned in the range +/- pi/2.
**
**  3)  The latitude phi is pi/2 minus the angle between the Earth's
**      rotation axis and the adopted zenith.  In many applications it
**      will be sufficient to use the published geodetic latitude of the
**      site.  In very precise (sub-arcsecond) applications, phi can be
**      corrected for polar motion.
**
**  4)  The returned azimuth az is with respect to the rotational north
**      pole, as opposed to the ITRS pole, and for sub-arcsecond
**      accuracy will need to be adjusted for polar motion if it is to
**      be with respect to north on a map of the Earth's surface.
**
**  5)  Should the user wish to work with respect to the astronomical
**      zenith rather than the geodetic zenith, phi will need to be
**      adjusted for deflection of the vertical (often tens of
**      arcseconds), and the zero point of the hour angle ha will also
**      be affected.
**
**  6)  The transformation is the same as Vh = Rz(pi)*Ry(pi/2-phi)*Ve,
**      where Vh and Ve are lefthanded unit vectors in the (az,el) and
**      (ha,dec) systems respectively and Ry and Rz are rotations about
**      first the y-axis and then the z-axis.  (n.b. Rz(pi) simply
**      reverses the signs of the x and y components.)  For efficiency,
**      the algorithm is written out rather than calling other utility
**      functions.  For applications that require even greater
**      efficiency, additional savings are possible if constant terms
**      such as functions of latitude are computed once and for all.
**
**  7)  Again for efficiency, no range checking of arguments is carried
**      out.
**
**  Last revision:   2017 September 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    az = c_double()
    el = c_double()
    _sofa.iauHd2ae(ha,dec,phi,byref(az),byref(el))
    return az.value, el.value

#iauHd2pa
try:
    _sofa.iauHd2pa.argtypes = [
    c_double, # ha
    c_double, # dec
    c_double, # phi
    ]
    _sofa.iauHd2pa.restype = c_double
except AttributeError:
    pass
def hd2pa(ha, dec, phi):
    """
    Parallactic angle for a given hour angle and declination.

    Parameters
    ----------
    ha : double
        hour angle.
    dec : double
        declination.
    phi : double
        parallactic angle.

    Returns
    -------
    double
        parallactic angle.
/*
**  - - - - - - - - -
**   i a u H d 2 p a
**  - - - - - - - - -
**
**  Parallactic angle for a given hour angle and declination.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     ha     double     hour angle
**     dec    double     declination
**     phi    double     site latitude
**
**  Returned (function value):
**            double     parallactic angle
**
**  Notes:
**
**  1)  All the arguments are angles in radians.
**
**  2)  The parallactic angle at a point in the sky is the position
**      angle of the vertical, i.e. the angle between the directions to
**      the north celestial pole and to the zenith respectively.
**
**  3)  The result is returned in the range -pi to +pi.
**
**  4)  At the pole itself a zero result is returned.
**
**  5)  The latitude phi is pi/2 minus the angle between the Earth's
**      rotation axis and the adopted zenith.  In many applications it
**      will be sufficient to use the published geodetic latitude of the
**      site.  In very precise (sub-arcsecond) applications, phi can be
**      corrected for polar motion.
**
**  6)  Should the user wish to work with respect to the astronomical
**      zenith rather than the geodetic zenith, phi will need to be
**      adjusted for deflection of the vertical (often tens of
**      arcseconds), and the zero point of the hour angle ha will also
**      be affected.
**
**  Reference:
**     Smart, W.M., "Spherical Astronomy", Cambridge University Press,
**     6th edition (Green, 1977), p49.
**
**  Last revision:   2017 September 12
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauHd2pa(ha,dec,phi)











##############################################################
#################################################################
########### Section2:    SOFA Vector/Matrix Library
######################################################
##############################################################   



#########################
### Initialize
############################
    
#iauZp
try:
    _sofa.iauZp.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3]
    ]
except AttributeError:
    pass
def zp():
    """
Zero a p-vector
    Returns
    -------
    p : double[3]
        p-vector.
/*
**  - - - - - -
**   i a u Z p
**  - - - - - -
**
**  Zero a p-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Returned:
**     p        double[3]      p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    p = zeros(shape=3,dtype = float)
    _sofa.iauZp(p)
    return  p

#iauZr
try:
    _sofa.iauZr.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ]
except AttributeError:
    pass
def zr():
    """
    
Initialize an r-matrix to the null matrix.
    Returns
    -------
    r : double[3][3]
         r-matrix.
/*
**  - - - - - -
**   i a u Z r
**  - - - - - -
**
**  Initialize an r-matrix to the null matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Returned:
**     r        double[3][3]    r-matrix
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = zeros(shape=(3,3),dtype = float)
    _sofa.iauZr(r)
    return  r


#iauir
try:
    _sofa.iauIr.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ]
except AttributeError:
    pass
def ir():
    """
    Initialize an r-matrix to the identity matrix.

    Returns
    -------
    r : double[3][3]
        r-matrix.
/*
**  - - - - - -
**   i a u I r
**  - - - - - -
**
**  Initialize an r-matrix to the identity matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Returned:
**     r       double[3][3]    r-matrix
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = zeros(shape=(3,3),dtype = float)
    _sofa.iauIr(r)
    return  r

##############################################
#####  Copy/extend/extract
#############################################

# iauCp   
try:
    _sofa.iauCp.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3][1]
    ndpointer(shape=3,dtype =float),#c[3][1]
    ]
except AttributeError:
    pass
def cp(p):
    """
    Copy a p-vector.

    Parameters
    ----------
    p : double[3]
        p-vector to be copied.

    Returns
    -------
    c : double[3]
        copy.
/*
**  - - - - - -
**   i a u C p
**  - - - - - -
**
**  Copy a p-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p        double[3]     p-vector to be copied
**
**  Returned:
**     c        double[3]     copy
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    c = zeros(shape=3,dtype = float)
    _sofa.iauCp(p,c)
    return  c

# iauCr   
try:
    _sofa.iauCr.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#p[3][3]
    ndpointer(shape=(3,3),dtype =float),#c[3][3]
    ]
except AttributeError:
    pass
def cr(p):
    """
    Copy an r-matrix.

    Parameters
    ----------
    p : double[3][3]
        r-matrix to be copied.

    Returns
    -------
    c : double[3][3]
        copy.
/*
**  - - - - - -
**   i a u C r
**  - - - - - -
**
**  Copy an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    r-matrix to be copied
**
**  Returned:
**     c        double[3][3]    copy
**
**  Called:
**     iauCp        copy p-vector
**
**  This revision:  2016 May 19
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    c = zeros(shape=(3,3),dtype = float)
    _sofa.iauCr(p,c)
    return  c

##################################################
####### Build rotations
##################################################
# iauRx  
try:
    _sofa.iauRx.argtypes = [
    c_double,#phi
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ]
except AttributeError:
    pass
def rx(phi,r):
    """
    Rotate an r-matrix about the x-axis.

    Parameters
    ----------
    phi : double
        angle (radians).
    r : doule[3][3]
        r-matrix, rotated.

    Returns
    -------
    r : TYPE
        DESCRIPTION.
/*
**  - - - - - -
**   i a u R x
**  - - - - - -
**
**  Rotate an r-matrix about the x-axis.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     phi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive phi incorporates in the
**     supplied r-matrix r an additional rotation, about the x-axis,
**     anticlockwise as seen looking towards the origin from positive x.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  1        0            0      )
**         (                               )
**         (  0   + cos(phi)   + sin(phi)  )
**         (                               )
**         (  0   - sin(phi)   + cos(phi)  )
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    _sofa.iauRx(phi,r)
    return  r

# iauRy 
try:
    _sofa.iauRy.argtypes = [
    c_double,#theta
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ]
except AttributeError:
    pass
def ry(theta,r):
    """
    Rotate an r-matrix about the y-axis.

    Parameters
    ----------
    theta : double
        angle (radians).
    r : double[3][3]
        r-matrix, rotated.

    Returns
    -------
    r : TYPE
        DESCRIPTION.
/*
**  - - - - - -
**   i a u R y
**  - - - - - -
**
**  Rotate an r-matrix about the y-axis.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     theta  double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive theta incorporates in the
**     supplied r-matrix r an additional rotation, about the y-axis,
**     anticlockwise as seen looking towards the origin from positive y.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(theta)     0      - sin(theta)  )
**         (                                        )
**         (       0           1           0        )
**         (                                        )
**         (  + sin(theta)     0      + cos(theta)  )
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    _sofa.iauRy(theta,r)
    return  r

# iauRz
try:
    _sofa.iauRz.argtypes = [
    c_double,#psi
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ]
except AttributeError:
    pass
def rz(psi,r):
    """
    Rotate an r-matrix about the z-axis.

    Parameters
    ----------
    psi : double
        angle (radians).
    r : double
        r-matrix, rotated.

    Returns
    -------
    r : double
        r-matrix, rotated.
/*
**  - - - - - -
**   i a u R z
**  - - - - - -
**
**  Rotate an r-matrix about the z-axis.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     psi    double          angle (radians)
**
**  Given and returned:
**     r      double[3][3]    r-matrix, rotated
**
**  Notes:
**
**  1) Calling this function with positive psi incorporates in the
**     supplied r-matrix r an additional rotation, about the z-axis,
**     anticlockwise as seen looking towards the origin from positive z.
**
**  2) The additional rotation can be represented by this matrix:
**
**         (  + cos(psi)   + sin(psi)     0  )
**         (                                 )
**         (  - sin(psi)   + cos(psi)     0  )
**         (                                 )
**         (       0            0         1  )
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    _sofa.iauRz(psi,r)
    return  r

######################################################
########  Spherical/Cartesian conversions
#######################################################
# iauS2c
try:
    _sofa.iauS2c.argtypes = [
    c_double,#theta
    c_double, # phi
    ndpointer(shape=3,dtype =float),#c[3]
    ]
except AttributeError:
    pass
def s2c(theta,phi):
    """
    Convert spherical coordinates to Cartesian.

    Parameters
    ----------
    theta : double
        longitude angle (radians).
    phi : double
        latitude angle (radians).

    Returns
    -------
    c : double[3]
        direction cosines.
/*
**  - - - - - - -
**   i a u S 2 c
**  - - - - - - -
**
**  Convert spherical coordinates to Cartesian.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**
**  Returned:
**     c        double[3]    direction cosines
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    c = zeros(shape=3,dtype=float)
    _sofa.iauS2c(theta,phi,c)
    return  c

# iauC2s
try:
    _sofa.iauC2s.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3]
    POINTER(c_double),#theta
    POINTER(c_double), # phi
    ]
except AttributeError:
    pass
def c2s(p):
    """
    P-vector to spherical coordinates.

    Parameters
    ----------
    p : double[3]
        p-vector.

    Returns
    -------
    theta : double
        longitude angle (radians).
    phi : double
        latitude angle (radians).
/*
**  - - - - - - -
**   i a u C 2 s
**  - - - - - - -
**
**  P-vector to spherical coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p      double[3]    p-vector
**
**  Returned:
**     theta  double       longitude angle (radians)
**     phi    double       latitude angle (radians)
**
**  Notes:
**
**  1) The vector p can have any magnitude; only its direction is used.
**
**  2) If p is null, zero theta and phi are returned.
**
**  3) At either pole, zero theta is returned.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    theta = c_double()
    phi = c_double()
    _sofa.iauC2s(p,byref(theta),byref(phi))
    return  theta.value, phi.value

# iauS2p
try:
    _sofa.iauS2p.argtypes = [
    c_double,#theta
    c_double, # phi
    c_double, # r
    ndpointer(shape=3,dtype =float),#p[3]
    ]
except AttributeError:
    pass
def s2p(theta,phi,r):
    """
    Convert spherical polar coordinates to p-vector.

    Parameters
    ----------
    theta : double
        longitude angle (radians).
    phi : double
        latitude angle (radians).
    r : double
        radial distance.

    Returns
    -------
    p : double[3]
        Cartesian coordinates.

    """
    p = zeros(shape=3,dtype=float)
    _sofa.iauS2p(theta,phi,r,p)
    return  p

# iauP2s
try:
    _sofa.iauP2s.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3]
    POINTER(c_double),#theta
    POINTER(c_double), # phi
    POINTER(c_double), #r
    ]
except AttributeError:
    pass
def p2s(p):
    """
    P-vector to spherical polar coordinates.

    Parameters
    ----------
    p : double[3]
        p-vector.

    Returns
    -------
    theta  : double
        longitude angle (radians).
    phi : double
        latitude angle (radians).
    r : double
        radial distance.
/*
**  - - - - - - -
**   i a u P 2 s
**  - - - - - - -
**
**  P-vector to spherical polar coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p        double[3]    p-vector
**
**  Returned:
**     theta    double       longitude angle (radians)
**     phi      double       latitude angle (radians)
**     r        double       radial distance
**
**  Notes:
**
**  1) If P is null, zero theta, phi and r are returned.
**
**  2) At either pole, zero theta is returned.
**
**  Called:
**     iauC2s       p-vector to spherical
**     iauPm        modulus of p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    theta = c_double()
    phi = c_double()
    r = c_double()
    _sofa.iauP2s(p,byref(theta),byref(phi),byref(r))
    return  theta.value, phi.value, r.value


#######################################################
##########  Operations on vectors
###################################################
# iauPpp
try:
    _sofa.iauPpp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ndpointer(shape=3,dtype =float),#apb[3]
    ]
except AttributeError:
    pass
def ppp(a,b):
    """
    P-vector addition.

    Parameters
    ----------
    a : double[3]
        first p-vector.
    b : double[3]
        second p-vector.

    Returns
    -------
    apb : double[3]
         a + b.
/*
**  - - - - - - -
**   i a u P p p
**  - - - - - - -
**
**  P-vector addition.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     apb      double[3]      a + b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    apb = zeros(shape=3,dtype = float)
    _sofa.iauPpp(a,b,apb)
    return  apb

# iauPmp
try:
    _sofa.iauPmp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ndpointer(shape=3,dtype =float),#amb[3]
    ]
except AttributeError:
    pass
def pmp(a,b):
    """
    P-vector subtraction.

    Parameters
    ----------
    a : double[3]
        first p-vector.
    b : double[3]
        second p-vector.

    Returns
    -------
    amb : TYPE
        DESCRIPTION.
/*
**  - - - - - - -
**   i a u P m p
**  - - - - - - -
**
**  P-vector subtraction.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     amb      double[3]      a - b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    amb = zeros(shape=3,dtype = float)
    _sofa.iauPmp(a,b,amb)
    return  amb

# iauPpsp
try:
    _sofa.iauPpsp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    c_double, # s
    ndpointer(shape=3,dtype =float),#b[3]
    ndpointer(shape=3,dtype =float),#apsb[3]
    ]
except AttributeError:
    pass
def ppsp(a,s,b):
    """
    P-vector plus scaled p-vector.

    Parameters
    ----------
    a : double[3]
        first p-vector.
    s : double
        scalar (multiplier for b).
    b : double[3]
        second p-vector.

    Returns
    -------
    apsb : TYPE
        DESCRIPTION.
/*
**  - - - - - - - -
**   i a u P p s p
**  - - - - - - - -
**
**  P-vector plus scaled p-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a      double[3]     first p-vector
**     s      double        scalar (multiplier for b)
**     b      double[3]     second p-vector
**
**  Returned:
**     apsb   double[3]     a + s*b
**
**  Note:
**     It is permissible for any of a, b and apsb to be the same array.
**
**  Called:
**     iauSxp       multiply p-vector by scalar
**     iauPpp       p-vector plus p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    apsb = zeros(shape=3,dtype = float)
    _sofa.iauPpsp(a,s,b,apsb)
    return  apsb

# iauPdp
try:
    _sofa.iauPdp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ]
    _sofa.iauPdp.restype = c_double
except AttributeError:
    pass
def pdp(a,b):
    """
    p-vector inner (=scalar=dot) product.

    Parameters
    ----------
    a : double[3]
        first p-vector.
    b : double[3]
        second p-vector.

    Returns
    -------
    ab : double
        a . b.
/*
**  - - - - - - -
**   i a u P d p
**  - - - - - - -
**
**  p-vector inner (=scalar=dot) product.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a      double[3]     first p-vector
**     b      double[3]     second p-vector
**
**  Returned (function value):
**            double        a . b
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    ab = _sofa.iauPdp(a,b)
    return  ab

# iauPxp
try:
    _sofa.iauPxp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ndpointer(shape=3,dtype =float),#axb[3]
    ]
except AttributeError:
    pass
def pxp(a,b):
    """
    p-vector outer (=vector=cross) product.

    Parameters
    ----------
    a : double[3]
        first p-vector.
    b : double[3]
        second p-vector.

    Returns
    -------
    axb : double[3]
         a x b.
/*
**  - - - - - - -
**   i a u P x p
**  - - - - - - -
**
**  p-vector outer (=vector=cross) product.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[3]      first p-vector
**     b        double[3]      second p-vector
**
**  Returned:
**     axb      double[3]      a x b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    axb = zeros(shape=3,dtype = float)
    _sofa.iauPxp(a,b,axb)
    return  axb

# iauPm
try:
    _sofa.iauPm.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3]
    ]
    _sofa.iauPm.restype = c_double
except AttributeError:
    pass
def pm(p):
    """
    Modulus of p-vector.

    Parameters
    ----------
    p : double[3]
        p-vector.

    Returns
    -------
    double
        modulus.
/*
**  - - - - - -
**   i a u P m
**  - - - - - -
**
**  Modulus of p-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p      double[3]     p-vector
**
**  Returned (function value):
**            double        modulus
**
**  This revision:  2013 August 7
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return  _sofa.iauPm(p)

# iauPn
try:
    _sofa.iauPn.argtypes = [
    ndpointer(shape=3,dtype =float),#p[3]
    POINTER(c_double),#r
    ndpointer(shape=3,dtype =float),#u[3]
    ]
except AttributeError:
    pass
def pn(p):
    """
    Convert a p-vector into modulus and unit vector.

    Parameters
    ----------
    p : double[3]
        p-vector.

    Returns
    -------
    r : double
        modulus.
    u : double[3]
        unit vector.
/*
**  - - - - - -
**   i a u P n
**  - - - - - -
**
**  Convert a p-vector into modulus and unit vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p        double[3]      p-vector
**
**  Returned:
**     r        double         modulus
**     u        double[3]      unit vector
**
**  Notes:
**
**  1) If p is null, the result is null.  Otherwise the result is a unit
**     vector.
**
**  2) It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     iauPm        modulus of p-vector
**     iauZp        zero p-vector
**     iauSxp       multiply p-vector by scalar
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = c_double()
    u = zeros(shape=3,dtype = float)
    _sofa.iauPn(p,byref(r),u)
    return  r.value, u

# iauSxp
try:
    _sofa.iauSxp.argtypes = [
    c_double,#s
    ndpointer(shape=3,dtype =float),#p[3]
    ndpointer(shape=3,dtype =float),#sp[3]
    ]
except AttributeError:
    pass
def sxp(s,p):
    """
    Multiply a p-vector by a scalar.

    Parameters
    ----------
    s : double
        scalar.
    p : double[3]
        p-vector.

    Returns
    -------
    sp : double[3]
        s*p.
/*
**  - - - - - - -
**   i a u S x p
**  - - - - - - -
**
**  Multiply a p-vector by a scalar.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     s      double        scalar
**     p      double[3]     p-vector
**
**  Returned:
**     sp     double[3]     s * p
**
**  Note:
**     It is permissible for p and sp to be the same array.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    sp = zeros(shape=3,dtype = float)
    _sofa.iauSxp(s,p,sp)
    return  sp


##############################################################
############### Operations on matrices
#####################################################
# iauRxr
try:
    _sofa.iauRxr.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#a[3][3]
    ndpointer(shape=(3,3),dtype =float),#b[3][3]
    ndpointer(shape=(3,3),dtype =float),#atb[3][3]
    ]
except AttributeError:
    pass
def rxr(a,b):
    """
    Multiply two r-matrices.

    Parameters
    ----------
    a : double[3][3]
        first r-matrix.
    b : double[3][3]
        second r-matrix.

    Returns
    -------
    atb : double[3][3]
        a * b.
/*
**  - - - - - - -
**   i a u R x r
**  - - - - - - -
**
**  Multiply two r-matrices.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[3][3]    first r-matrix
**     b        double[3][3]    second r-matrix
**
**  Returned:
**     atb      double[3][3]    a * b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     iauCr        copy r-matrix
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    atb = zeros(shape=(3,3),dtype = float)
    _sofa.iauRxr(a,b,atb)
    return  atb

# iauTr
try:
    _sofa.iauTr.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ndpointer(shape=(3,3),dtype =float),#b[3][3]
    ]
except AttributeError:
    pass
def tr(r):
    """
    Transpose an r-matrix.

    Parameters
    ----------
    r : double[3][3]
        r-matrix.

    Returns
    -------
    rt : double[3][3]
        transpose.
/*
**  - - - - - -
**   i a u T r
**  - - - - - -
**
**  Transpose an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    r-matrix
**
**  Returned:
**     rt       double[3][3]    transpose
**
**  Note:
**     It is permissible for r and rt to be the same array.
**
**  Called:
**     iauCr        copy r-matrix
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rt = zeros(shape=(3,3),dtype = float)
    _sofa.iauTr(r,rt)
    return  rt


##############################################################
############### Matrix−vector products
#####################################################
# iauRxp
try:
    _sofa.iauRxp.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ndpointer(shape=3,dtype =float),#p[3]
    ndpointer(shape=3,dtype =float),#rp[3]
    ]
except AttributeError:
    pass
def rxp(r,p):
    """
    Multiply a p-vector by an r-matrix.

    Parameters
    ----------
    r : double[3][3]
        r-matrix.
    p : double[3]
        p-vector.

    Returns
    -------
    rp : double[3]
        r * p.
/*
**  - - - - - - -
**   i a u R x p
**  - - - - - - -
**
**  Multiply a p-vector by an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    r-matrix
**     p        double[3]       p-vector
**
**  Returned:
**     rp       double[3]       r * p
**
**  Note:
**     It is permissible for p and rp to be the same array.
**
**  Called:
**     iauCp        copy p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rp = zeros(shape=3,dtype = float)
    _sofa.iauRxp(r,p,rp)
    return rp

# iauTrxp
try:
    _sofa.iauTrxp.argtypes = [
    ndpointer(shape=(3,3),dtype =float),#r[3][3]
    ndpointer(shape=3,dtype =float),#p[3]
    ndpointer(shape=3,dtype =float),#trp[3]
    ]
except AttributeError:
    pass
def trxp(r,p):
    """
    Multiply a p-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : double[3][3]
        r-matrix.
    p : double[3]
        p-vector.

    Returns
    -------
    trp : double[3]
        r * p.
/*
**  - - - - - - - -
**   i a u T r x p
**  - - - - - - - -
**
**  Multiply a p-vector by the transpose of an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]   r-matrix
**     p        double[3]      p-vector
**
**  Returned:
**     trp      double[3]      r * p
**
**  Note:
**     It is permissible for p and trp to be the same array.
**
**  Called:
**     iauTr        transpose r-matrix
**     iauRxp       product of r-matrix and p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    trp = zeros(shape=3,dtype = float)
    _sofa.iauTrxp(r,p,trp)
    return trp


##############################################################
############### Separation and position−angle 
#####################################################
# iauSepp
try:
    _sofa.iauSepp.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ]
    _sofa.iauSepp.restype = c_double
except AttributeError:
    pass
def sepp(a,b):
    """
    Angular separation between two p-vectors.

    Parameters
    ----------
    a : double[3]
        first p-vector (not necessarily unit length).
    b : double[3]
        second p-vector (not necessarily unit length).

    Returns
    -------
    double
        angular separation (radians, always positive).
/*
**  - - - - - - - -
**   i a u S e p p
**  - - - - - - - -
**
**  Angular separation between two p-vectors.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a      double[3]    first p-vector (not necessarily unit length)
**     b      double[3]    second p-vector (not necessarily unit length)
**
**  Returned (function value):
**            double       angular separation (radians, always positive)
**
**  Notes:
**
**  1) If either vector is null, a zero result is returned.
**
**  2) The angular separation is most simply formulated in terms of
**     scalar product.  However, this gives poor accuracy for angles
**     near zero and pi.  The present algorithm uses both cross product
**     and dot product, to deliver full accuracy whatever the size of
**     the angle.
**
**  Called:
**     iauPxp       vector product of two p-vectors
**     iauPm        modulus of p-vector
**     iauPdp       scalar product of two p-vectors
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauSepp(a,b)



# iauSeps
try:
    _sofa.iauSeps.argtypes = [
    c_double,#al
    c_double,#ap
    c_double,#bl
    c_double,#bp
    ]
    _sofa.iauSeps.restype = c_double
except AttributeError:
    pass
def seps(al,ap,bl,bp):
    """
    Angular separation between two sets of spherical coordinates.

    Parameters
    ----------
    al : double
        first longitude (radians).
    ap : double
        first latitude (radians).
    bl : double
        second longitude (radians).
    bp : double
        second latitude (radians).

    Returns
    -------
    double
        angular separation (radians).
/*
**  - - - - - - - -
**   i a u S e p s
**  - - - - - - - -
**
**  Angular separation between two sets of spherical coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     al     double       first longitude (radians)
**     ap     double       first latitude (radians)
**     bl     double       second longitude (radians)
**     bp     double       second latitude (radians)
**
**  Returned (function value):
**            double       angular separation (radians)
**
**  Called:
**     iauS2c       spherical coordinates to unit vector
**     iauSepp      angular separation between two p-vectors
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauSeps(al,ap,bl,bp)


# iauPap
try:
    _sofa.iauPap.argtypes = [
    ndpointer(shape=3,dtype =float),#a[3]
    ndpointer(shape=3,dtype =float),#b[3]
    ]
    _sofa.iauPap.restype = c_double
except AttributeError:
    pass
def pap(a,b):
    """
    Position-angle from two p-vectors.

    Parameters
    ----------
    a : double[3]
        direction of reference point.
    b : double[3]
        direction of point whose PA is required.

    Returns
    -------
    double
        position angle of b with respect to a (radians).
/*
**  - - - - - - -
**   i a u P a p
**  - - - - - - -
**
**  Position-angle from two p-vectors.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a      double[3]  direction of reference point
**     b      double[3]  direction of point whose PA is required
**
**  Returned (function value):
**            double     position angle of b with respect to a (radians)
**
**  Notes:
**
**  1) The result is the position angle, in radians, of direction b with
**     respect to direction a.  It is in the range -pi to +pi.  The
**     sense is such that if b is a small distance "north" of a the
**     position angle is approximately zero, and if b is a small
**     distance "east" of a the position angle is approximately +pi/2.
**
**  2) The vectors a and b need not be of unit length.
**
**  3) Zero is returned if the two directions are the same or if either
**     vector is null.
**
**  4) If vector a is at a pole, the result is ill-defined.
**
**  Called:
**     iauPn        decompose p-vector into modulus and direction
**     iauPm        modulus of p-vector
**     iauPxp       vector product of two p-vectors
**     iauPmp       p-vector minus p-vector
**     iauPdp       scalar product of two p-vectors
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauPap(a,b)


# iauPas
try:
    _sofa.iauPas.argtypes = [
    c_double,#al
    c_double,#ap
    c_double,#bl
    c_double,#bp
    ]
    _sofa.iauPas.restype = c_double
except AttributeError:
    pass
def pas(al,ap,bl,bp):
    """
    Position-angle from spherical coordinates.

    Parameters
    ----------
    al : double
        longitude of point A (e.g. RA) in radians.
    ap : double
        latitude of point A (e.g. Dec) in radians.
    bl : double
        longitude of point B.
    bp : double
        latitude of point B.

    Returns
    -------
    double
        position angle of B with respect to A.
/*
**  - - - - - - -
**   i a u P a s
**  - - - - - - -
**
**  Position-angle from spherical coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     al     double     longitude of point A (e.g. RA) in radians
**     ap     double     latitude of point A (e.g. Dec) in radians
**     bl     double     longitude of point B
**     bp     double     latitude of point B
**
**  Returned (function value):
**            double     position angle of B with respect to A
**
**  Notes:
**
**  1) The result is the bearing (position angle), in radians, of point
**     B with respect to point A.  It is in the range -pi to +pi.  The
**     sense is such that if B is a small distance "east" of point A,
**     the bearing is approximately +pi/2.
**
**  2) Zero is returned if the two points are coincident.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauPas(al,ap,bl,bp)


##############################################################
############### Rotation vectors
#####################################################
 
#iauRv2m
try:
    _sofa.iauRv2m.argtypes = [
    ndpointer(shape=3,dtype = float),#w[3]
    ndpointer(shape=(3,3),dtype = float),#r[3][3]
    ]
except AttributeError:
    pass
def rv2m(w):
    """
    
Form the r-matrix corresponding to a given r-vector.
    Parameters
    ----------
    w : double[3]
        rotation vector (Note 1).

    Returns
    -------
    r : double[3][3]
        rotation matrix.
/*
**  - - - - - - - -
**   i a u R v 2 m
**  - - - - - - - -
**
**  Form the r-matrix corresponding to a given r-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     w        double[3]      rotation vector (Note 1)
**
**  Returned:
**     r        double[3][3]    rotation matrix
**
**  Notes:
**
**  1) A rotation matrix describes a rotation through some angle about
**     some arbitrary axis called the Euler axis.  The "rotation vector"
**     supplied to This function has the same direction as the Euler
**     axis, and its magnitude is the angle in radians.
**
**  2) If w is null, the unit matrix is returned.
**
**  3) The reference frame rotates clockwise as seen looking along the
**     rotation vector from the origin.
**
**  This revision:  2015 January 30
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = zeros(shape=(3,3),dtype = float)
    _sofa.iauRv2m(w,r)
    return r

#iauRm2v
try:
    _sofa.iauRm2v.argtypes = [
    ndpointer(shape=(3,3),dtype = float),#r[3][3]
    ndpointer(shape=3,dtype = float),#w[3][1]
    ]
except AttributeError:
    pass
def rm2v(r):
    """
    Express an r-matrix as an r-vector.

    Parameters
    ----------
    r : double[3][3]
        rotation matrix.

    Returns
    -------
    w : double[3]
        rotation vector (Note 1).
/*
**  - - - - - - - -
**   i a u R m 2 v
**  - - - - - - - -
**
**  Express an r-matrix as an r-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    rotation matrix
**
**  Returned:
**     w        double[3]       rotation vector (Note 1)
**
**  Notes:
**
**  1) A rotation matrix describes a rotation through some angle about
**     some arbitrary axis called the Euler axis.  The "rotation vector"
**     returned by this function has the same direction as the Euler axis,
**     and its magnitude is the angle in radians.  (The magnitude and
**     direction can be separated by means of the function iauPn.)
**
**  2) If r is null, so is the result.  If r is not a rotation matrix
**     the result is undefined;  r must be proper (i.e. have a positive
**     determinant) and real orthogonal (inverse = transpose).
**
**  3) The reference frame rotates clockwise as seen looking along
**     the rotation vector from the origin.
**
**  This revision:  2015 January 30
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    w = zeros(shape=3,dtype = float)
    _sofa.iauRm2v(r,w)
    return w



######################################
#########################################################################
############### Section 2.2 OPERATIONS INVOLVING PV−VECTORS
#########################################################################
#####################
    
##############################################################
###############  Initialize
#####################################################  
#iauZpv
try:
    _sofa.iauZpv.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ]
except AttributeError:
    pass
def zpv():
    """
    Zero a pv-vector.

    Returns
    -------
    pv : double[2][3]
         pv-vector.
/*
**  - - - - - - -
**   i a u Z p v
**  - - - - - - -
**
**  Zero a pv-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Returned:
**     pv       double[2][3]      pv-vector
**
**  Called:
**     iauZp        zero p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    pv = zeros(shape=(2,3),dtype = float)
    _sofa.iauZpv(pv)
    return pv
##############################################################
###############  Copy/extend/extract
#####################################################  
#iauCpv
try:
    _sofa.iauCpv.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#c[2][3]
    ]
except AttributeError:
    pass
def cpv(pv):
    """
    Copy a position/velocity vector.

    Parameters
    ----------
    pv : double[2][3]
        position/velocity vector to be copied.

    Returns
    -------
    c : double[2][3]
        copy.
/*
**  - - - - - - -
**   i a u C p v
**  - - - - - - -
**
**  Copy a position/velocity vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     pv     double[2][3]    position/velocity vector to be copied
**
**  Returned:
**     c      double[2][3]    copy
**
**  Called:
**     iauCp        copy p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    c = zeros(shape=(2,3),dtype = float)
    _sofa.iauCpv(pv,c)
    return c

#iauP2pv
try:
    _sofa.iauP2pv.argtypes = [
    ndpointer(shape=3,dtype = float),#p[3]
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ]
except AttributeError:
    pass
def p2pv(p):
    """
    Extend a p-vector to a pv-vector by appending a zero velocity.

    Parameters
    ----------
    p : doube[3]
        p-vector.

    Returns
    -------
    pv : double[2][3]
        pv-vector.
/*
**  - - - - - - - -
**   i a u P 2 p v
**  - - - - - - - -
**
**  Extend a p-vector to a pv-vector by appending a zero velocity.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     p        double[3]       p-vector
**
**  Returned:
**     pv       double[2][3]    pv-vector
**
**  Called:
**     iauCp        copy p-vector
**     iauZp        zero p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    pv = zeros(shape=(2,3),dtype = float)
    _sofa.iauP2pv(p,pv)
    return pv

#iauPv2p
try:
    _sofa.iauPv2p.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#pv[3]
    ndpointer(shape=3,dtype = float),#p[3]
    ]
except AttributeError:
    pass
def pv2p(pv):
    """
    Discard velocity component of a pv-vector.

    Parameters
    ----------
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    p : double[3]
        p-vector.
/*
**  - - - - - - - -
**   i a u P v 2 p
**  - - - - - - - -
**
**  Discard velocity component of a pv-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     pv      double[2][3]     pv-vector
**
**  Returned:
**     p       double[3]        p-vector
**
**  Called:
**     iauCp        copy p-vector
**
**  This revision:  2013 June 18	
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    p = zeros(shape=3,dtype = float)
    _sofa.iauPv2p(pv,p)
    return p


##############################################################
###############  Spherical/Cartesian conversions
##################################################### 
#iauS2pv
try:
    _sofa.iauS2pv.argtypes = [
    c_double,#theta
    c_double,# phi
    c_double, # r
    c_double, #td
    c_double, #pd
    c_double, # rd
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ]
except AttributeError:
    pass
def s2pv(theta,phi,r,td,pd,rd):
    """
    Convert position/velocity from spherical to Cartesian coordinates.

    Parameters
    ----------
    theta : double
        longitude angle (radians).
    phi : double
        latitude angle (radians).
    r : double
        radial distance.
    td : double
        rate of change of theta.
    pd : double
        rate of change of phi.
    rd : double
        rate of change of r.

    Returns
    -------
    pv : double[2][3]
        pv-vector.
/*
**  - - - - - - - -
**   i a u S 2 p v
**  - - - - - - - -
**
**  Convert position/velocity from spherical to Cartesian coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     theta    double          longitude angle (radians)
**     phi      double          latitude angle (radians)
**     r        double          radial distance
**     td       double          rate of change of theta
**     pd       double          rate of change of phi
**     rd       double          rate of change of r
**
**  Returned:
**     pv       double[2][3]    pv-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    pv = zeros(shape=(2,3),dtype = float)
    _sofa.iauS2pv(theta,phi,r,td,pd,rd,pv)
    return pv

#iauPv2s
try:
    _sofa.iauPv2s.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    POINTER(c_double),#theta
    POINTER(c_double),# phi
    POINTER(c_double), # r
    POINTER(c_double), #td
    POINTER(c_double), #pd
    POINTER(c_double), # rd
    
    ]
except AttributeError:
    pass
def pv2s(pv):
    """
    Convert position/velocity from Cartesian to spherical coordinates.

    Parameters
    ----------
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    theta : double
        longitude angle (radians).
    phi : double
         latitude angle (radians).
    r : double
        radial distance.
    td : double
        rate of change of theta.
    pd : double
        rate of change of phi.
    rd : double
        rate of change of r.
/*
**  - - - - - - - -
**   i a u P v 2 s
**  - - - - - - - -
**
**  Convert position/velocity from Cartesian to spherical coordinates.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     pv       double[2][3]  pv-vector
**
**  Returned:
**     theta    double        longitude angle (radians)
**     phi      double        latitude angle (radians)
**     r        double        radial distance
**     td       double        rate of change of theta
**     pd       double        rate of change of phi
**     rd       double        rate of change of r
**
**  Notes:
**
**  1) If the position part of pv is null, theta, phi, td and pd
**     are indeterminate.  This is handled by extrapolating the
**     position through unit time by using the velocity part of
**     pv.  This moves the origin without changing the direction
**     of the velocity component.  If the position and velocity
**     components of pv are both null, zeroes are returned for all
**     six results.
**
**  2) If the position is a pole, theta, td and pd are indeterminate.
**     In such cases zeroes are returned for all three.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    theta = c_double()
    phi = c_double()
    r = c_double()
    td  = c_double()
    pd = c_double()
    rd = c_double()
    _sofa.iauPv2s(pv,byref(theta),byref(phi),byref(r),byref(td),byref(pd),byref(rd))
    return theta.value, phi.value, r.value,td.value, pd.value, rd.value


##############################################################
###############  Operations on vectors
#####################################################  
#iauPvppv
try:
    _sofa.iauPvppv.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#a[2][3]
    ndpointer(shape=(2,3),dtype = float),#b[2][3]
    ndpointer(shape=(2,3),dtype = float),#apb[2][3]
    ]
except AttributeError:
    pass
def pvppv(a,b):
    """
    Add one pv-vector to another.

    Parameters
    ----------
    a : double[2][3]
        first pv-vector.
    b : double[2][3]
        second pv-vector.

    Returns
    -------
    apb : double[2][3]
        a+b.
/*
**  - - - - - - - - -
**   i a u P v p p v
**  - - - - - - - - -
**
**  Add one pv-vector to another.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[2][3]      first pv-vector
**     b        double[2][3]      second pv-vector
**
**  Returned:
**     apb      double[2][3]      a + b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     iauPpp       p-vector plus p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    apb = zeros(shape=(2,3), dtype = float)
    _sofa.iauPvppv(a,b,apb)
    return apb

#iauPvmpv
try:
    _sofa.iauPvmpv.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#a[2][3]
    ndpointer(shape=(2,3),dtype = float),#b[2][3]
    ndpointer(shape=(2,3),dtype = float),#amb[2][3]
    ]
except AttributeError:
    pass
def pvmpv(a,b):
    """
    Subtract one pv-vector from another.

    Parameters
    ----------
    a : double[2][3]
        first pv-vector.
    b : dobule[2][3]
        second pv-vector.

    Returns
    -------
    amb : double[2][3]
        a - b.
/*
**  - - - - - - - - -
**   i a u P v m p v
**  - - - - - - - - -
**
**  Subtract one pv-vector from another.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a       double[2][3]      first pv-vector
**     b       double[2][3]      second pv-vector
**
**  Returned:
**     amb     double[2][3]      a - b
**
**  Note:
**     It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     iauPmp       p-vector minus p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    amb = zeros(shape=(2,3), dtype = float)
    _sofa.iauPvmpv(a,b,amb)
    return amb

#iauPvxpv
try:
    _sofa.iauPvxpv.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#a[2][3]
    ndpointer(shape=(2,3),dtype = float),#b[2][3]
    ndpointer(shape=(2,3),dtype = float),#axb[2][3]
    ]
except AttributeError:
    pass
def pvxpv(a,b):
    """
    Outer (=vector=cross) product of two pv-vectors.

    Parameters
    ----------
    a : double[2][3]
        first pv-vector.
    b : double[2][3]
        second pv-vector.

    Returns
    -------
    axb : double[2][3]
        a *b.
/*
**  - - - - - - - - -
**   i a u P v x p v
**  - - - - - - - - -
**
**  Outer (=vector=cross) product of two pv-vectors.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double[2][3]      first pv-vector
**     b        double[2][3]      second pv-vector
**
**  Returned:
**     axb      double[2][3]      a x b
**
**  Notes:
**
**  1) If the position and velocity components of the two pv-vectors are
**     ( ap, av ) and ( bp, bv ), the result, a x b, is the pair of
**     vectors ( ap x bp, ap x bv + av x bp ).  The two vectors are the
**     cross-product of the two p-vectors and its derivative.
**
**  2) It is permissible to re-use the same array for any of the
**     arguments.
**
**  Called:
**     iauCpv       copy pv-vector
**     iauPxp       vector product of two p-vectors
**     iauPpp       p-vector plus p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    axb = zeros(shape=(2,3), dtype = float)
    _sofa.iauPvxpv(a,b,axb)
    return axb

#iauPvm
try:
    _sofa.iauPvm.argtypes = [
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    POINTER(c_double), #r
    POINTER(c_double), #s
    ]
except AttributeError:
    pass
def pvm(pv):
    """
    Modulus of pv-vector.

    Parameters
    ----------
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    r : double
        modulus of position component.
    s : double
        modulus of velocity component.
/*
**  - - - - - - -
**   i a u P v m
**  - - - - - - -
**
**  Modulus of pv-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     pv     double[2][3]   pv-vector
**
**  Returned:
**     r      double         modulus of position component
**     s      double         modulus of velocity component
**
**  Called:
**     iauPm        modulus of p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    r = c_double()
    s = c_double()
    _sofa.iauPvm(pv,byref(r),byref(s))
    return r.value, s.value

#iauSxpv
try:
    _sofa.iauSxpv.argtypes = [
    c_double,#s
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#spv[2][3]
    ]
except AttributeError:
    pass
def sxpv(s,pv):
    """
    Multiply a pv-vector by a scalar.

    Parameters
    ----------
    s : double
        scalar.
    pv : double
        pv-vector.

    Returns
    -------
    spv : double[2][3]
        s * pv.
/*
**  - - - - - - - -
**   i a u S x p v
**  - - - - - - - -
**
**  Multiply a pv-vector by a scalar.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     s       double          scalar
**     pv      double[2][3]    pv-vector
**
**  Returned:
**     spv     double[2][3]    s * pv
**
**  Note:
**     It is permissible for pv and spv to be the same array
**
**  Called:
**     iauS2xpv     multiply pv-vector by two scalars
**
**  This revision:  2013 August 7
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    spv = zeros(shape=(2,3), dtype = float)
    _sofa.iauSxpv(s,pv,spv)
    return spv

#iauS2xpv
try:
    _sofa.iauS2xpv.argtypes = [
    c_double,#s1
    c_double,#s2
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#spv[2][3]
    ]
except AttributeError:
    pass
def s2xpv(s1,s2,pv):
    """
    Multiply a pv-vector by two scalars.

    Parameters
    ----------
    s1 : double
        scalar to multiply position component by.
    s2 : double
        scalar to multiply velocity component by.
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    spv : double[2][3]
        pv-vector: p scaled by s1, v scaled by s2.
/*
**  - - - - - - - - -
**   i a u S 2 x p v
**  - - - - - - - - -
**
**  Multiply a pv-vector by two scalars.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     s1     double         scalar to multiply position component by
**     s2     double         scalar to multiply velocity component by
**     pv     double[2][3]   pv-vector
**
**  Returned:
**     spv    double[2][3]   pv-vector: p scaled by s1, v scaled by s2
**
**  Note:
**     It is permissible for pv and spv to be the same array.
**
**  Called:
**     iauSxp       multiply p-vector by scalar
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    spv = zeros(shape=(2,3), dtype = float)
    _sofa.iauS2xpv(s1,s2,pv,spv)
    return spv

#iauPvu
try:
    _sofa.iauPvu.argtypes = [
    c_double,#dt
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#upv[2][3]
    ]
except AttributeError:
    pass
def pvu(dt,pv):
    """
    Update a pv-vector.

    Parameters
    ----------
    dt : double
        time interval.
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    upv : double[2][3]
        p updated, v unchanged.
/*
**  - - - - - - -
**   i a u P v u
**  - - - - - - -
**
**  Update a pv-vector.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     dt       double           time interval
**     pv       double[2][3]     pv-vector
**
**  Returned:
**     upv      double[2][3]     p updated, v unchanged
**
**  Notes:
**
**  1) "Update" means "refer the position component of the vector
**     to a new date dt time units from the existing date".
**
**  2) The time units of dt must match those of the velocity.
**
**  3) It is permissible for pv and upv to be the same array.
**
**  Called:
**     iauPpsp      p-vector plus scaled p-vector
**     iauCp        copy p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    upv = zeros(shape=(2,3), dtype = float)
    _sofa.iauPvu(dt,pv,upv)
    return upv

#iauPvup
try:
    _sofa.iauPvup.argtypes = [
    c_double,#dt
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=3,dtype = float),#p[3]
    ]
except AttributeError:
    pass
def pvup(dt,pv):
    """
    Update a pv-vector, discarding the velocity component.

    Parameters
    ----------
    dt : double
        time interval.
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    p : double[3]
        p-vector.
/*
**  - - - - - - - -
**   i a u P v u p
**  - - - - - - - -
**
**  Update a pv-vector, discarding the velocity component.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     dt       double            time interval
**     pv       double[2][3]      pv-vector
**
**  Returned:
**     p        double[3]         p-vector
**
**  Notes:
**
**  1) "Update" means "refer the position component of the vector to a
**     new date dt time units from the existing date".
**
**  2) The time units of dt must match those of the velocity.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    p = zeros(shape=3, dtype = float)
    _sofa.iauPvup(dt,pv,p)
    return p



##############################################################
###############  Matrix−vector products
#####################################################  
#iauRxpv
try:
    _sofa.iauRxpv.argtypes = [
    ndpointer(shape=(3,3),dtype = float),#r[3][3]
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#rpv[2][3]
    ]
    
except AttributeError:
    pass
def rxpv(r,pv):
    """
    Multiply a pv-vector by an r-matrix.

    Parameters
    ----------
    r : double[3][3]
        r-matrix.
    pv : double[2][3]
        pv-vector.

    Returns
    -------
    rpv : double[2][3]
        r * pv.
/*
**  - - - - - - - -
**   i a u R x p v
**  - - - - - - - -
**
**  Multiply a pv-vector by an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    r-matrix
**     pv       double[2][3]    pv-vector
**
**  Returned:
**     rpv      double[2][3]    r * pv
**
**  Note:
**     It is permissible for pv and rpv to be the same array.
**
**  Called:
**     iauRxp       product of r-matrix and p-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    rpv = zeros(shape=(2,3),dtype = float)
    _sofa.iauRxpv(r,pv,rpv)
    return rpv

#iauTrxpv
try:
    _sofa.iauTrxpv.argtypes = [
    ndpointer(shape=(3,3),dtype = float),#r[3][3]
    ndpointer(shape=(2,3),dtype = float),#pv[2][3]
    ndpointer(shape=(2,3),dtype = float),#trpv[2][3]
    ]
    
except AttributeError:
    pass
def trxpv(r,pv):
    """
    Multiply a pv-vector by the transpose of an r-matrix.

    Parameters
    ----------
    r : double[3][3]
        r-matrix.
    pv : double[3][3]
        pv-vector.

    Returns
    -------
    trpv : double[2][3]
         r * pv.
/*
**  - - - - - - - - -
**   i a u T r x p v
**  - - - - - - - - -
**
**  Multiply a pv-vector by the transpose of an r-matrix.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     r        double[3][3]    r-matrix
**     pv       double[2][3]    pv-vector
**
**  Returned:
**     trpv     double[2][3]    r * pv
**
**  Note:
**     It is permissible for pv and trpv to be the same array.
**
**  Called:
**     iauTr        transpose r-matrix
**     iauRxpv      product of r-matrix and pv-vector
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    trpv = zeros(shape=(2,3),dtype = float)
    _sofa.iauTrxpv(r,pv,trpv)
    return trpv
    

##############################################################
###############   OPERATIONS ON ANGLES
#####################################################  

#iauAnp
try:
    _sofa.iauAnp.argtypes = [
    c_double, # a
    ]
    _sofa.iauAnp.restype = c_double
    
except AttributeError:
    pass
def anp(a):
    """
    Normalize angle into the range 0 <= a < 2pi.

    Parameters
    ----------
    a : double
        angle (radians).

    Returns
    -------
    double
        angle in range 0-2pi.
/*
**  - - - - - - -
**   i a u A n p
**  - - - - - - -
**
**  Normalize angle into the range 0 <= a < 2pi.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range 0-2pi
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauAnp(a)


#iauAnpm
try:
    _sofa.iauAnpm.argtypes = [
    c_double, # a
    ]
    _sofa.iauAnpm.restype = c_double
    
except AttributeError:
    pass
def anpm(a):
    """
    Normalize angle into the range -pi <= a < +pi.

    Parameters
    ----------
    a : double
        angle(radians).

    Returns
    -------
    double
        angle in range +/-pi.
/*
**  - - - - - - - -
**   i a u A n p m
**  - - - - - - - -
**
**  Normalize angle into the range -pi <= a < +pi.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     a        double     angle (radians)
**
**  Returned (function value):
**              double     angle in range +/-pi
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    return _sofa.iauAnpm(a)

#iauA2tf
try:
    _sofa.iauA2tf.argtypes = [
    c_int, # ndp
    c_double, # angle
    POINTER(c_char), #sign
    ndpointer(shape = (4,),dtype = int), #ihmsf
    ]
    
except AttributeError:
    pass
def a2tf(ndp, angle):
    """
    Decompose radians into hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution (Note 1).
    angle : double
        angle in radians.

    Returns
    -------
    sign : char
        '+' or '-'.
    ihmsf : TYPE
        hours, minutes, seconds, fraction.
/*
**  - - - - - - - -
**   i a u A 2 t f
**  - - - - - - - -
**
**  Decompose radians into hours, minutes, seconds, fraction.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     angle   double  angle in radians
**
**  Returned:
**     sign    char    '+' or '-'
**     ihmsf   int[4]  hours, minutes, seconds, fraction
**
**  Called:
**     iauD2tf      decompose days to hms
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of angle, the format of doubles on the target platform, and
**     the risk of overflowing ihmsf[3].  On a typical platform, for
**     angle up to 2pi, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of angle may exceed 2pi.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where angle is very nearly 2pi and rounds up to 24 hours,
**     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
**
**  This revision:  2013 July 31
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    sign = c_char()
    ihmsf = zeros(shape = (4,),dtype = int)
    _sofa.iauA2tf(ndp,angle,byref(sign),ihmsf)
    return sign.value.decode(), ihmsf

#iauA2af
try:
    _sofa.iauA2af.argtypes = [
    c_int, # ndp
    c_double, # angle
    POINTER(c_char), #sign
    ndpointer(shape = (4,),dtype = int), #idmsf
    ]
    
except AttributeError:
    pass
def a2af(ndp, angle):
    """
    Decompose radians into degrees, arcminutes, arcseconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution (Note 1).
    angle : double
        angle in radians.

    Returns
    -------
    sign : char
        '+' or '-'.
    idmsf : int[4]
        degrees, arcminutes, arcseconds, fraction.
/*
**  - - - - - - - -
**   i a u A 2 a f
**  - - - - - - - -
**
**  Decompose radians into degrees, arcminutes, arcseconds, fraction.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     angle   double  angle in radians
**
**  Returned:
**     sign    char    '+' or '-'
**     idmsf   int[4]  degrees, arcminutes, arcseconds, fraction
**
**  Called:
**     iauD2tf      decompose days to hms
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of angle, the format of doubles on the target platform, and
**     the risk of overflowing idmsf[3].  On a typical platform, for
**     angle up to 2pi, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of angle may exceed 2pi.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where angle is very nearly 2pi and rounds up to 360 degrees,
**     by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    sign = c_char()
    idmsf = zeros(shape = (4,),dtype = int)
    _sofa.iauA2af(ndp,angle,byref(sign),idmsf)
    return sign.value.decode(), idmsf

#iauAf2a
try:
    _sofa.iauAf2a.argtypes = [
    c_char, # s
    c_int, # ideg
    c_int, # iamin
    c_double, # asec
    POINTER(c_double), #rad
    ]
    
except AttributeError:
    pass
def af2a(s, ideg, iamin, asec):
    """
    Convert degrees, arcminutes, arcseconds to radians.

    Parameters
    ----------
    s : char
        sign:  '-' = negative, otherwise positive.
    ideg : int
        degrees.
    iamin : int
        arcminutes.
    asec : double
        arcseconds.

    Returns
    -------
    rad : double
        double  angle in radians.
/*
**  - - - - - - - -
**   i a u A f 2 a
**  - - - - - - - -
**
**  Convert degrees, arcminutes, arcseconds to radians.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     s         char    sign:  '-' = negative, otherwise positive
**     ideg      int     degrees
**     iamin     int     arcminutes
**     asec      double  arcseconds
**
**  Returned:
**     rad       double  angle in radians
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = ideg outside range 0-359
**                                2 = iamin outside range 0-59
**                                3 = asec outside range 0-59.999...
**
**  Notes:
**
**  1)  The result is computed even if any of the range checks fail.
**
**  2)  Negative ideg, iamin and/or asec produce a warning status, but
**      the absolute value is used in the conversion.
**
**  3)  If there are multiple errors, the status value reflects only the
**      first, the smallest taking precedence.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = s.encode()
    rad = c_double()
    _sofa.iauAf2a(s, ideg, iamin, asec, byref(rad))
    return rad.value

#iauD2tf
try:
    _sofa.iauD2tf.argtypes = [
    c_int, # ndp
    c_double, # days
    POINTER(c_char), #sign
    ndpointer(shape = (4,),dtype = int), #ihmsf
    ]
    
except AttributeError:
    pass
def d2tf(ndp, days):
    """
    Decompose days to hours, minutes, seconds, fraction.

    Parameters
    ----------
    ndp : int
        resolution (Note 1).
    days : double
        interval in days.

    Returns
    -------
    sign : char
        '+' or '-'.
    ihmsf : TYPE
        hours, minutes, seconds, fraction.
/*
**  - - - - - - - -
**   i a u D 2 t f
**  - - - - - - - -
**
**  Decompose days to hours, minutes, seconds, fraction.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards Of Fundamental Astronomy) software collection.
**
**  Status:  vector/matrix support function.
**
**  Given:
**     ndp     int     resolution (Note 1)
**     days    double  interval in days
**
**  Returned:
**     sign    char    '+' or '-'
**     ihmsf   int[4]  hours, minutes, seconds, fraction
**
**  Notes:
**
**  1) The argument ndp is interpreted as follows:
**
**     ndp         resolution
**      :      ...0000 00 00
**     -7         1000 00 00
**     -6          100 00 00
**     -5           10 00 00
**     -4            1 00 00
**     -3            0 10 00
**     -2            0 01 00
**     -1            0 00 10
**      0            0 00 01
**      1            0 00 00.1
**      2            0 00 00.01
**      3            0 00 00.001
**      :            0 00 00.000...
**
**  2) The largest positive useful value for ndp is determined by the
**     size of days, the format of double on the target platform, and
**     the risk of overflowing ihmsf[3].  On a typical platform, for
**     days up to 1.0, the available floating-point precision might
**     correspond to ndp=12.  However, the practical limit is typically
**     ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
**     only 16 bits.
**
**  3) The absolute value of days may exceed 1.0.  In cases where it
**     does not, it is up to the caller to test for and handle the
**     case where days is very nearly 1.0 and rounds up to 24 hours,
**     by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    sign = c_char()
    ihmsf = zeros(shape = (4,),dtype = int)
    _sofa.iauD2tf(ndp,days,byref(sign),ihmsf)
    return sign.value.decode(), ihmsf

#iauTf2a
try:
    _sofa.iauTf2a.argtypes = [
    c_char, # s
    c_int, # ihour
    c_int, # imin
    c_double, # sec
    POINTER(c_double), #rad
    ]
    
except AttributeError:
    pass
def tf2a(s, ihour, imin, sec):
    """
    Convert hours, minutes, seconds to radians.

    Parameters
    ----------
    s : char
        sign:  '-' = negative, otherwise positive.
    ihour : int
        hours.
    imin : int
        minutes.
    sec : double
        seconds.

    Returns
    -------
    rad : double
        angle in radians.
/*
**  - - - - - - - -
**   i a u T f 2 a
**  - - - - - - - -
**
**  Convert hours, minutes, seconds to radians.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     s         char    sign:  '-' = negative, otherwise positive
**     ihour     int     hours
**     imin      int     minutes
**     sec       double  seconds
**
**  Returned:
**     rad       double  angle in radians
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = ihour outside range 0-23
**                                2 = imin outside range 0-59
**                                3 = sec outside range 0-59.999...
**
**  Notes:
**
**  1)  The result is computed even if any of the range checks fail.
**
**  2)  Negative ihour, imin and/or sec produce a warning status, but
**      the absolute value is used in the conversion.
**
**  3)  If there are multiple errors, the status value reflects only the
**      first, the smallest taking precedence.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    s = s.encode()
    rad = c_double()
    _sofa.iauTf2a(s, ihour, imin, sec,rad)
    return rad.value

#iauTf2d
try:
    _sofa.iauTf2d.argtypes = [
    c_char, # s
    c_int, # ihour
    c_int, # imin
    c_double, # sec
    POINTER(c_double), #days
    ]
    
except AttributeError:
    pass
def tf2d(s, ihour, imin, sec):
    """
    Convert hours, minutes, seconds to days.

    Parameters
    ----------
    s : char
        sign:  '-' = negative, otherwise positive.
    ihour : int
        hours.
    imin : int
        minutes.
    sec : double
        seconds.

    Returns
    -------
    days : double
        interval in days.
/*
**  - - - - - - - -
**   i a u T f 2 d
**  - - - - - - - -
**
**  Convert hours, minutes, seconds to days.
**
**  This function is part of the International Astronomical Union's
**  SOFA (Standards of Fundamental Astronomy) software collection.
**
**  Status:  support function.
**
**  Given:
**     s         char    sign:  '-' = negative, otherwise positive
**     ihour     int     hours
**     imin      int     minutes
**     sec       double  seconds
**
**  Returned:
**     days      double  interval in days
**
**  Returned (function value):
**               int     status:  0 = OK
**                                1 = ihour outside range 0-23
**                                2 = imin outside range 0-59
**                                3 = sec outside range 0-59.999...
**
**  Notes:
**
**  1)  The result is computed even if any of the range checks fail.
**
**  2)  Negative ihour, imin and/or sec produce a warning status, but
**      the absolute value is used in the conversion.
**
**  3)  If there are multiple errors, the status value reflects only the
**      first, the smallest taking precedence.
**
**  This revision:  2013 June 18
**
**  SOFA release 2019-07-22
**
**  Copyright (C) 2019 IAU SOFA Board.  See notes at end.
*/
    """
    days = c_double()
    s = s.encode()
    _sofa.iauTf2d(s, ihour, imin, sec,days)
    return days.value
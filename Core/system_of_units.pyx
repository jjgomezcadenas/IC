# $Id: system_of_units.h,v 1.1 2004/04/16 17:09:03 gomez Exp $
# ----------------------------------------------------------------------
# HEP coherent system of Unitscoulomb
#
# This file has been provided by Geant4 (simulation toolkit for HEP).
#
# The basic units are :
#  		millimeter              (millimeter)
# 		nanosecond              (nanosecond)
# 		Mega electron Volt      (MeV)
# 		positron charge         (eplus)
# 		degree Kelvin           (kelvin)
#              the amount of substance (mole)
#              luminous intensity      (candela)
# 		radian                  (radian)
#              steradian               (steradian)
#
# Below is a non exhaustive list of derived and pratical units
# (i.e. mostly the SI units).
# You can add your own units.
#
# The SI numerical value of the positron charge is defined here,
# as it is needed for conversion factor : positron charge = e_SI (coulomb)
#
# The others physical constants are defined in the header file :
# PhysicalConstants.h
#
# Authors: M.Maire, S.Giani
#
# History:
#
# 06.02.96   Created.
# 28.03.96   Added miscellaneous constants.
# 05.12.97   E.Tcherniaev: Redefined pascal (to avoid warnings on WinNT)
# 20.05.98   names: meter, second, gram, radian, degree
#            (from Brian.Lasiuk@yale.edu (STAR)). Added luminous units.
# 05.08.98   angstrom, picobarn, microsecond, picosecond, petaelectronvolt


#
# Length [L]
#

cdef class SystemOfUnits:
    cdef float euro = 1.
    cdef float millimeter = 1.
    cdef float millimeter2 = millimeter * millimeter
    cdef float millimeter3 = millimeter * millimeter2

    cdef float centimeter = 10. * millimeter
    cdef float centimeter2 = centimeter * centimeter
    cdef float centimeter3 = centimeter * centimeter2

    cdef float decimeter = 100. * millimeter
    cdef float decimeter2 = decimeter * decimeter
    cdef float decimeter3 = decimeter * decimeter2
    cdef float liter = decimeter3
    cdef float l = liter

    cdef float meter = 1000. * millimeter
    cdef float meter2 = meter * meter
    cdef float meter3 = meter * meter2

    cdef float kilometer = 1000. * meter
    cdef float kilometer2 = kilometer * kilometer
    cdef float kilometer3 = kilometer * kilometer2

    cdef float micrometer = 1.e-6 * meter
    cdef float nanometer = 1.e-9 * meter
    cdef float angstrom = 1.e-10 * meter
    cdef float fermi = 1.e-15 * meter

    cdef float nm = nanometer
    cdef float mum = micrometer

    cdef float micron = micrometer
    cdef float barn = 1.e-28 * meter2
    cdef float millibarn = 1.e-3 * barn
    cdef float microbarn = 1.e-6 * barn
    cdef float nanobarn = 1.e-9 * barn
    cdef float picobarn = 1.e-12 * barn

    # symbols
    cdef float mm = millimeter
    cdef float mm2 = millimeter2
    cdef float mm3 = millimeter3

    cdef float cm = centimeter
    cdef float cm2 = centimeter2
    cdef float cm3 = centimeter3

    cdef float m = meter
    cdef float m2 = meter2
    cdef float m3 = meter3

    cdef float km = kilometer
    cdef float km2 = kilometer2
    cdef float km3 = kilometer3

    cdef float ft = 30.48 * cm

    #
    # Angle
    #
    cdef float radian = 1.
    cdef float milliradian = 1.e-3 * radian
    cdef float degree = (3.14159265358979323846/180.0) * radian

    cdef float steradian = 1.

    # symbols
    cdef float rad = radian
    cdef float mrad = milliradian
    cdef float sr = steradian
    cdef float deg = degree

    #
    # Time [T]
    #
    cdef float nanosecond = 1.
    cdef float second = 1.e+9 * nanosecond
    cdef float millisecond = 1.e-3 * second
    cdef float microsecond = 1.e-6 * second
    cdef float picosecond = 1.e-12 * second
    cdef float year = 3.1536e+7 * second
    cdef float day = 864e2 * second
    cdef float minute = 60 * second
    cdef float hour = 60 * minute

    cdef float s = second
    cdef float ms = millisecond
    cdef float ps = picosecond
    cdef float mus = microsecond
    cdef float ns = nanosecond

    # new!
    cdef float hertz = 1./second
    cdef float kilohertz = 1.e+3 * hertz
    cdef float megahertz = 1.e+6 * hertz
    cdef float gigahertz = 1.e+6 * hertz

    cdef float MHZ = megahertz
    cdef float kHZ = kilohertz
    cdef float kHz = kHZ
    cdef float GHZ = gigahertz

    #
    # Electric charge [Q]
    #
    cdef float eplus = 1.  # positron charge
    cdef float e_SI = 1.60217733e-19  # positron charge in coulomb
    cdef float coulomb = eplus/e_SI  # coulomb = 6.24150 e+18 * eplus
    #
    # Energy [E]
    #
    cdef float megaelectronvolt = 1.
    cdef float electronvolt = 1.e-6 * megaelectronvolt
    cdef float milielectronvolt = 1.e-3 * electronvolt
    cdef float kiloelectronvolt = 1.e-3 * megaelectronvolt
    cdef float gigaelectronvolt = 1.e+3 * megaelectronvolt
    cdef float teraelectronvolt = 1.e+6 * megaelectronvolt
    cdef float petaelectronvolt = 1.e+9 * megaelectronvolt

    cdef float meV = milielectronvolt
    cdef float eV = electronvolt
    cdef float keV = kiloelectronvolt
    cdef float MeV = megaelectronvolt
    cdef float GeV = gigaelectronvolt
    cdef float TeV = teraelectronvolt
    cdef float PeV = petaelectronvolt
    cdef float eV2 = eV*eV

    cdef float joule = electronvolt/e_SI  # joule = 6.24150 e+12 * MeV
    #
    # Mass [E][T^2][L^-2]
    #
    cdef float kilogram = joule * second * second / meter2
    cdef float gram = 1.e-3 * kilogram
    cdef float milligram = 1.e-3 * gram
    cdef float ton = 1.e+3 * kilogram
    cdef float kiloton = 1.e+3 * ton

    # symbols
    cdef float kg = kilogram
    cdef float g = gram
    cdef float mg = milligram
    #
    # Power [E][T^-1]
    #
    cdef float watt = joule/second  # watt = 6.24150 e+3 * MeV/ns
    #
    # Force [E][L^-1]
    #
    cdef float newton = joule/meter  # newton = 6.24150 e+9 * MeV/mm
    #
    # Pressure [E][L^-3]
    #
    cdef float hep_pascal = newton/m2  # pascal = 6.24150 e+3 * MeV/mm3
    cdef float pascal = hep_pascal
    cdef float Pa = pascal
    cdef float kPa = 1000 * Pa
    cdef float MPa = 1e+6 * Pa
    cdef float GPa = 1e+9 * Pa
    cdef float bar = 100000 * pascal  # bar = 6.24150 e+8 * MeV/mm3
    cdef float milibar = 1e-3 * bar
    cdef float atmosphere = 101325 * pascal  # atm = 6.32420 e+8 * MeV/mm3
    cdef float denier = gram / (9000 * meter)
    #
    # Electric current [Q][T^-1]
    cdef float ampere = coulomb/second  # ampere = 6.24150 e+9 * eplus/ns
    cdef float milliampere = 1.e-3 * ampere
    cdef float microampere = 1.e-6 * ampere
    cdef float nanoampere = 1.e-9 * ampere
    cdef float mA = milliampere
    cdef float muA = microampere
    cdef float nA = nanoampere
    #
    # Electric potential [E][Q^-1]
    #
    cdef float megavolt = megaelectronvolt/eplus
    cdef float kilovolt = 1.e-3 * megavolt
    cdef float volt = 1.e-6 * megavolt
    cdef float millivolt = 1.e-3 * volt

    cdef float V = volt
    cdef float mV = millivolt
    cdef float kV = kilovolt
    cdef float MV = megavolt

    #
    # Electric resistance [E][T][Q^-2]
    #
    cdef float ohm = volt/ampere  # ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)
    #
    # Electric capacitance [Q^2][E^-1]
    #
    cdef float farad = coulomb/volt  # farad = 6.24150e+24 * eplus/Megavolt
    cdef float millifarad = 1.e-3 * farad
    cdef float microfarad = 1.e-6 * farad
    cdef float nanofarad = 1.e-9 * farad
    cdef float picofarad = 1.e-12 * farad

    cdef float nF = nanofarad
    cdef float pF = picofarad
    #
    # Magnetic Flux [T][E][Q^-1]
    #
    cdef float weber = volt * second  # weber = 1000*megavolt*ns
    #
    # Magnetic Field [T][E][Q^-1][L^-2]
    #
    cdef float tesla = volt*second/meter2  # tesla = 0.001*megavolt*ns/mm2
    cdef float gauss = 1.e-4 * tesla
    cdef float kilogauss = 1.e-1 * tesla
    #
    # Inductance [T^2][E][Q^-2]
    #
    cdef float henry = weber/ampere  # henry = 1.60217e-7*MeV*(ns/eplus)**2
    #
    # Temperature
    #
    cdef float kelvin = 1.
    #
    # Amount of substance
    #
    cdef float mole = 1.
    cdef float mol = mole
    #
    # Activity [T^-1]
    #
    cdef float becquerel = 1./second
    cdef float curie = 3.7e+10 * becquerel
    cdef float Bq = becquerel
    cdef float mBq = 1e-3 * becquerel
    cdef float muBq = 1e-6 * becquerel
    cdef float cks = Bq/keV
    cdef float U238ppb = Bq / 81.
    cdef float Th232ppb = Bq / 246.
    #
    # Absorbed dose [L^2][T^-2]
    #
    cdef float gray = joule/kilogram
    #
    # Luminous intensity [I]
    #
    cdef float candela = 1.
    #
    # Luminous flux [I]
    #
    cdef float lumen = candela * steradian
    #
    # Illuminance [I][L^-2]
    #
    cdef float lux = lumen/meter2
    #
    # Miscellaneous
    #
    cdef float perCent = 1e-2
    cdef float perThousand = 1e-3
    cdef float perMillion = 1e-6
    cdef float pes = 1.
    cdef float adc = 1

cpdef celsius(float tKelvin):
    return tKelvin - 273.15

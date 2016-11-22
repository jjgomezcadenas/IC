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
    cdef readonly double euro, millimeter, millimeter2, millimeter3, centimeter, centimeter2, cm, cm2, cm3, centimeter3, decimeter, decimeter2, decimeter3, liter, l, meter, meter2, meter3, kilometer,  kilometer2,  kilometer3,  micrometer,  nanometer,  angstrom, fermi,  nm,  mum,  micron,  barn,  millibarn,  microbarn,  nanobarn,  picobarn,  mm,  mm2,  mm3, m,  m2,  m3,  km, km2, km3, ft,  radian,  milliradian, mrad, degree,  steradian,  rad, sr,  deg,  nanosecond,  millisecond,  second, year,  day,  minute,  hour,  s,  ms,  ps, mus, ns, picosecond, microsecond, hertz, kilohertz, megahertz,  gigahertz,  MHZ,  kHz, kHZ, GHZ,  eplus,  e_SI,  coulomb,  electronvolt, megaelectronvolt,  milielectronvolt,  kiloelectronvolt, gigaelectronvolt,  teraelectronvolt, petaelectronvolt, meV,  eV,  keV,  MeV,  GeV,  TeV,  PeV, eV2,  joule,  kilogram,  gram,  milligram, ton,  kiloton,  kg,  g,  mg,  watt, newton, hep_pascal,  pascal,  Pa,  kPa,  MPa,  GPa,  bar,  milibar,  atmosphere,  denier,  ampere,  milliampere, microampere, nanoampere,  mA,  muA,  nA,  megavolt,  kilovolt,  volt,  millivolt,  V,  mV,  kV, MV,  ohm, farad,  millifarad,  microfarad,  nanofarad,  picofarad, nF,  pF, weber,  tesla,  gauss,  kilogauss,  henry,  kelvin,  mole,  mol,  becquerel,  curie,  Bq,  mBq,  muBq, cks,  U238ppb, Th232ppb,  gray,  candela,  lumen,  lux, perCent, perThousand,  perMillion,  pes, adc


    def __init__(self):
        self.euro = 1.
        self.millimeter = 1.
        self.millimeter2 = self.millimeter * self.millimeter
        self.millimeter3 = self.millimeter * self.millimeter2
#
        self.centimeter = 10. * self.millimeter
        self.centimeter2 = self.centimeter * self.centimeter
        self.centimeter3 = self.centimeter * self.centimeter2
#
        self.decimeter = 100. * self.millimeter
        self.decimeter2 = self.decimeter * self.decimeter
        self.decimeter3 = self.decimeter * self.decimeter2
        self.liter = self.decimeter3
        self.l = self.liter
#
        self.meter = 1000. * self.millimeter
        self.meter2 = self.meter * self.meter
        self.meter3 = self.meter * self.meter2
#
        self.kilometer = 1000. * self.meter
        self.kilometer2 = self.kilometer * self.kilometer
        self.kilometer3 = self.kilometer * self.kilometer2
#
        self.micrometer = 1.e-6 * self.meter
        self.nanometer = 1.e-9 * self.meter
        self.angstrom = 1.e-10 * self.meter
        self.fermi = 1.e-15 * self.meter
#
        self.nm = self.nanometer
        self.mum = self.micrometer
#
        self.micron = self.micrometer
        self.barn = 1.e-28 * self.meter2
        self.millibarn = 1.e-3 * self.barn
        self.microbarn = 1.e-6 * self.barn
        self.nanobarn = 1.e-9 * self.barn
        self.picobarn = 1.e-12 * self.barn
#
#     # symbols
        self.mm = self.millimeter
        self.mm2 = self.millimeter2
        self.mm3 = self.millimeter3
#
        self.cm = self.centimeter
        self.cm2 = self.centimeter2
        self.cm3 = self.centimeter3
#
        self.m = self.meter
        self.m2 = self.meter2
        self.m3 = self.meter3
#
        self.km = self.kilometer
        self.km2 = self.kilometer2
        self.km3 = self.kilometer3
#
        self.ft = 30.48 * self.cm
#
#     #
#     # Angle
#     #
        self.radian = 1.
        self.milliradian = 1.e-3 * self.radian
        self.degree = (3.14159265358979323846/180.0) * self.radian
#
        self.steradian = 1.
#
#     # symbols
        self.rad = self.radian
        self.mrad = self.milliradian
        self.sr = self.steradian
        self.deg = self.degree
#
#     #
#     # Time [T]
#     #
        self.nanosecond = 1.
        self.second = 1.e+9 * self.nanosecond
        self.millisecond = 1.e-3 * self.second
        self.microsecond = 1.e-6 * self.second
        self.picosecond = 1.e-12 * self.second
        self.year = 3.1536e+7 * self.second
        self.day = 864e2 * self.second
        self.minute = 60 * self.second
        self.hour = 60 * self.minute
#
        self.s = self.second
        self.ms = self.millisecond
        self.ps = self.picosecond
        self.mus = self.microsecond
        self.ns = self.nanosecond
#
#     # new!
        self.hertz = 1./self.second
        self.kilohertz = 1.e+3 * self.hertz
        self.megahertz = 1.e+6 * self.hertz
        self.gigahertz = 1.e+6 * self.hertz
#
        self.MHZ = self.megahertz
        self.kHZ = self.kilohertz
        self.kHz = self.kHZ
        self.GHZ = self.gigahertz
#
#     #
#     # Electric charge [Q]
#     #
        self.eplus = 1.  # positron charge
        self.e_SI = 1.60217733e-19  # positron charge in coulomb
        self.coulomb = self.eplus/self.e_SI  # coulomb = 6.24150 e+18 * eplus
#     #
#     # Energy [E]
#     #
        self.megaelectronvolt = 1.
        self.electronvolt = 1.e-6 * self.megaelectronvolt
        self.milielectronvolt = 1.e-3 * self.electronvolt
        self.kiloelectronvolt = 1.e-3 * self.megaelectronvolt
        self.gigaelectronvolt = 1.e+3 * self.megaelectronvolt
        self.teraelectronvolt = 1.e+6 * self.megaelectronvolt
        self.petaelectronvolt = 1.e+9 * self.megaelectronvolt
#
        self.meV = self.milielectronvolt
        self.eV = self.electronvolt
        self.keV = self.kiloelectronvolt
        self.MeV = self.megaelectronvolt
        self.GeV = self.gigaelectronvolt
        self.TeV = self.teraelectronvolt
        self.PeV = self.petaelectronvolt
        self.eV2 = self.eV*self.eV
#
        self.joule = self.electronvolt/self.e_SI  # joule = 6.24150 e+12 * MeV
#     #
#     # Mass [E][T^2][L^-2]
#     #
        self.kilogram = self.joule * self.second * self.second / self.meter2
        self.gram = 1.e-3 * self.kilogram
        self.milligram = 1.e-3 * self.gram
        self.ton = 1.e+3 * self.kilogram
        self.kiloton = 1.e+3 * self.ton
#
#     # symbols
        self.kg = self.kilogram
        self.g = self.gram
        self.mg = self.milligram
#     #
#     # Power [E][T^-1]
#     #
        self.watt = self.joule/self.second  # watt = 6.24150 e+3 * MeV/ns
#     #
#     # Force [E][L^-1]
#     #
        self.newton = self.joule/self.meter  # newton = 6.24150 e+9 * MeV/mm
#     #
#     # Pressure [E][L^-3]
#     #
        self.hep_pascal = self.newton/self.m2  # pascal = 6.24150 e+3 * MeV/mm3
        self.pascal = self.hep_pascal
        self.Pa = self.pascal
        self.kPa = 1000 * self.Pa
        self.MPa = 1e+6 * self.Pa
        self.GPa = 1e+9 * self.Pa
        self.bar = 100000 * self.pascal  # bar = 6.24150 e+8 * MeV/mm3
        self.milibar = 1e-3 * self.bar
        self.atmosphere = 101325 * self.pascal  # atm = 6.32420 e+8 * MeV/mm3
        self.denier = self.gram / (9000 * self.meter)
#     #
#     # Electric current [Q][T^-1]
        self.ampere = self.coulomb/self.second  # ampere = 6.24150 e+9 * eplus/ns
        self.milliampere = 1.e-3 * self.ampere
        self.microampere = 1.e-6 * self.ampere
        self.nanoampere = 1.e-9 * self.ampere
        self.mA = self.milliampere
        self.muA = self.microampere
        self.nA = self.nanoampere
#     #
#     # Electric potential [E][Q^-1]
#     #
        self.megavolt = self.megaelectronvolt/self.eplus
        self.kilovolt = 1.e-3 * self.megavolt
        self.volt = 1.e-6 * self.megavolt
        self.millivolt = 1.e-3 * self.volt
#
        self.V = self.volt
        self.mV = self.millivolt
        self.kV = self.kilovolt
        self.MV = self.megavolt
#
#     #
#     # Electric resistance [E][T][Q^-2]
#     #
        self.ohm = self.volt/self.ampere  # ohm = 1.60217e-16*(MeV/eplus)/(eplus/ns)
#     #
#     # Electric capacitance [Q^2][E^-1]
#     #
        self.farad = self.coulomb/self.volt  # farad = 6.24150e+24 * eplus/Megavolt
        self.millifarad = 1.e-3 * self.farad
        self.microfarad = 1.e-6 * self.farad
        self.nanofarad = 1.e-9 * self.farad
        self.picofarad = 1.e-12 * self.farad
#
        self.nF = self.nanofarad
        self.pF = self.picofarad
#     #
#     # Magnetic Flux [T][E][Q^-1]
#     #
        self.weber = self.volt * self.second  # weber = 1000*megavolt*ns
#     #
#     # Magnetic Field [T][E][Q^-1][L^-2]
#     #
        self.tesla = self.volt*self.second/self.meter2  # tesla = 0.001*megavolt*ns/mm2
        self.gauss = 1.e-4 * self.tesla
        self.kilogauss = 1.e-1 * self.tesla
#     #
#     # Inductance [T^2][E][Q^-2]
#     #
        self.henry = self.weber/self.ampere  # henry = 1.60217e-7*MeV*(ns/eplus)**2
#     #
#     # Temperature
#     #
        self.kelvin = 1.
#     #
#     # Amount of substance
#     #
        self.mole = 1.
        self.mol = self.mole
#     #
#     # Activity [T^-1]
#     #
        self.becquerel = 1./self.second
        self.curie = 3.7e+10 * self.becquerel
        self.Bq = self.becquerel
        self.mBq = 1e-3 * self.becquerel
        self.muBq = 1e-6 * self.becquerel
        self.cks = self.Bq/self.keV
        self.U238ppb = self.Bq / 81.
        self.Th232ppb = self.Bq / 246.
#     #
#     # Absorbed dose [L^2][T^-2]
#     #
        self.gray = self.joule/self.kilogram
#     #
#     # Luminous intensity [I]
#     #
        self.candela = 1.
#     #
#     # Luminous flux [I]
#     #
        self.lumen = self.candela * self.steradian
#     #
#     # Illuminance [I][L^-2]
#     #
        self.lux = self.lumen/self.meter2
#     #
#     # Miscellaneous
#     #
        self.perCent = 1e-2
        self.perThousand = 1e-3
        self.perMillion = 1e-6
        self.pes = 1.
        self.adc = 1
#
cpdef celsius(float tKelvin):
    return tKelvin - 273.15

from __future__ import division
import core

system = core.system('S.I.')

m = core.base_unit('length', 'meter', 'm', system)
kg = core.base_unit('mass', 'kilogram', 'kg', system)
s = core.base_unit('time', 'second', 's' , system)
A = core.base_unit('electric current', 'ampere', 'A', system)
K = core.base_unit('thermodynamic temperature ', 'kelvin', 'K', system)
mol = core.base_unit('amount of substance', 'mole', 'mol', system)
cd = core.base_unit('luminous intensity', 'candela', 'cd', system)

no_unit = core.none_unit('no unit','()')
none = core.none_unit('no unit','()')

named = core.named_unit

#SI prefixes
def deca(unit):
    """si prefix"""
    return core.factor_unit(10.,unit,name='deca'+unit.name,symbol='da'+unit.symbol)    
def hecto(unit):
    """si prefix"""
    return core.factor_unit(100.,unit,name='hecto'+unit.name,symbol='h'+unit.symbol)    
def kilo(unit):
    """si prefix"""
    return core.factor_unit(1000.,unit,name='kilo'+unit.name,symbol='k'+unit.symbol)    
def mega(unit):
    """si prefix"""
    return core.factor_unit(1.e6,unit,name='mega'+unit.name,symbol='M'+unit.symbol)    
def giga(unit):
    """si prefix"""
    return core.factor_unit(1.e9,unit,name='giga'+unit.name,symbol='G'+unit.symbol)    
def terra(unit):
    """si prefix"""
    return core.factor_unit(1.e12,unit,name='terra'+unit.name,symbol='T'+unit.symbol)    
def peta(unit):
    """si prefix"""
    return core.factor_unit(1.e15,unit,name='peta'+unit.name,symbol='P'+unit.symbol)    
def exa(unit):
    """si prefix"""
    return core.factor_unit(1.e18,unit,name='exa'+unit.name,symbol='E'+unit.symbol)    
def zetta(unit):
    """si prefix"""
    return core.factor_unit(1.e21,unit,name='zetta'+unit.name,symbol='Z'+unit.symbol)    
def yotta(unit):
    """si prefix"""
    return core.factor_unit(1.e24,unit,name='yotta'+unit.name,symbol='Y'+unit.symbol)    
def deci(unit):
    """si prefix"""
    return core.factor_unit(0.1,unit,name='deci'+unit.name,symbol='d'+unit.symbol)    
def centi(unit):
    """si prefix"""
    return core.factor_unit(0.01,unit,name='centi'+unit.name,symbol='c'+unit.symbol)    
def milli(unit):
    """si prefix"""
    return core.factor_unit(0.001,unit,name='milli'+unit.name,symbol='m'+unit.symbol)    
def micro(unit):
    """si prefix"""
    return core.factor_unit(1.e-6,unit,name='micro'+unit.name,symbol='mu'+unit.symbol)    
def nano(unit):
    """si prefix"""
    return core.factor_unit(1.e-9,unit,name='nano'+unit.name,symbol='n'+unit.symbol)    
def pico(unit):
    """si prefix"""
    return core.factor_unit(1.e-12,unit,name='pico'+unit.name,symbol='p'+unit.symbol)    
def femto(unit):
    """si prefix"""
    return core.factor_unit(1.e-15,unit,name='femto'+unit.name,symbol='f'+unit.symbol)    
def atto(unit):
    """si prefix"""
    return core.factor_unit(1.e-18,unit,name='atto'+unit.name,symbol='a'+unit.symbol)    
def zepto(unit):
    """si prefix"""
    return core.factor_unit(1.e-21,unit,name='zepto'+unit.name,symbol='z'+unit.symbol)    
def yocto(unit):
    """si prefix"""
    return core.factor_unit(1.e-24,unit,name='yocto'+unit.name,symbol='y'+unit.symbol)    
    


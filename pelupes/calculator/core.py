from __future__ import division
import values

class late(object):    
    """
    An attribute that is set at first access. 
    
    """ 
    def __init__(self, initializer):
        
        self.initializer = initializer
        self.__doc__ = self.initializer.__doc__
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        value = self.initializer(instance)
        setattr(instance,self.initializer.__name__, value)
        return value

class system(object):
    def __init__(self, name):
        self.name = name
        self.bases = []
    
    def add_base(self, unit):
        unit.system = self
        unit.index = len(self.bases)
        self.bases.append(unit)
        
class unit(object):
    """
    Abstract base class for unit objects.
        
    Two classes of units are defined:

    base units
        The base units in a given system of units. For SI, these
        are meter, kilogram, second, ampere, kelvin, mole and 
        candele. See the si module :mod:`amuse.support.units.si`
    derived units
        Derived units are created by dividing or multiplying
        with a number or with another unit. For example, 
        to get a velocity unit we can devine vel = 1000 * m / s

    Units can also be named, by creating a named unit.
    """
    
    def __mul__(self, other):
        if isinstance(other, unit):
#            if isinstance(other.base_unit(),none_unit):
#              return other.factor*self
#            if isinstance(self.base_unit(),none_unit):
#              return self.factor*other
            return mul_unit(self, other)
        else:
            return other*self
#            return factor_unit(other, self)
        
    def __div__(self, other):
        if isinstance(other, unit):
#            if isinstance(other.base_unit(),none_unit):
#              return (1./other.factor)*self
#            if isinstance(self.base_unit(),none_unit):
#              return self.factor/other
#            return self * other**-1
            return div_unit(self, other)
        else:
            return (1.0/other)*self
#            return factor_unit(1.0 / other, self)
    __truediv__=__div__

    def __rmul__(self, other):
        if other == 1.0:
            return self
        else:
            return factor_unit(other, self)

    def __or__(self, other):
        if isinstance(other, unit):
          return mul_unit(self,other)
        if isinstance(other, values.Quantity):
          return (1|self)*other
        return self.new_quantity(other)
    
    def __ror__(self, other):
        """Create a new Quantity object.
    
        :argument value: numeric value of the quantity, can be 
            a number or a sequence (list or ndarray)
        :returns: new ScalarQuantity or VectorQuantity object 
            with this unit
            
        Examples
        
        >>> from amuse.support.units import units
        >>> 100 | units.kg
        quantity<100 kg>
        """
        if isinstance(other, unit):
          return mul_unit(self,other)
        if isinstance(other, values.Quantity):
          return other*(1|self)
        return self.new_quantity(other)

    @property
    def number(self):
        return 1.0
        
    @property
    def unit(self):
        return self
        
    def __rdiv__(self, other):
#        return factor_unit(other, pow_unit(-1,self))
        return other * self**-1
        
    def __pow__(self, other):
        if other == 1:
            return self
        else:
            return pow_unit(other, self)
        
    def __call__(self, x):
        return self.new_quantity(x)
        
    
    def new_quantity(self, value):
        """Create a new Quantity object.
    
        :argument value: numeric value of the quantity, can be 
            a number or a sequence (list or ndarray)
        :returns: new ScalarQuantity or VectorQuantity object 
            with this unit
        """
        return values.ScalarQuantity(value,self) 
        
    def to_simple_form(self):
        """Convert unit to a form with only one factor and powers
        
        :result: Unit with only a factor and power terms
        
        >>> from amuse.support.units import units
        >>> N = (units.m * units.kg) / (units.s * units.s)
        >>> N
        unit<m * kg / s * s>
        >>> J = N * units.m
        >>> J
        unit<m * kg / s * s * m>
        >>> J.to_simple_form()
        unit<m**2 * kg * s**-2>
        """
        
        if not self.base:
            return none_unit('none', 'none') * self.factor
        
        result = 1
        for n, base in self.base:
            if n == 1:
                if result == 1:
                    result = base
                else:
                    result =  result*base 
            else:
                result = result*(base ** n) 
        result = self.factor*result
        
        return result
    
    def are_bases_equal(self, other):
        for n1, unit1 in sort(self.base, key=lambda x: x[1].index):
            found = False
            for n2, unit2 in other.base:
                if unit1 == unit2:
                    if not n2 == n1:
                        return False
                    found = True
                    break;
            if not found:
                return False
        return True
                        
    def conversion_factor_from(self, x):
        if self.base == x.base:
            this_factor = self.factor * 1.0
            other_factor = x.factor
            return this_factor / other_factor
        else:
            raise Exception("Cannot expres: " + str(x) + " in " + str(self))
      
    def in_(self, x):
        return self.as_quantity_in(x)

    def convert_to(self, x):
        return self.as_quantity_in(x)
    
    def as_quantity_in(self, unit):
        """Express this unit as a quantity in the given unit
        
        :argument unit: The unit to express this unit in
        :result: A Quantity object
        
        Examples
        
        >>> from amuse.support.units import units
        >>> ton = 1000 * units.kg
        >>> ton.as_quantity_in(units.kg)
        quantity<1000.0 kg>
        """
        if isinstance(unit, values.Quantity):
            raise Exception("Cannot expres a unit in a quantity")
        else:
            factor = self.conversion_factor_from(unit)
            return factor | unit
            
    def value_in(self, unit):
        """
        Return a numeric value of this unit in the given unit.
        Works only when the units are compatible, i.e. from
        tonnage to kg's.
        
        A number is returned without any unit information.
        
        :argument unit: wanted unit of the value
        :returns: number in the given unit
        
        >>> from amuse.support.units import units
        >>> x = units.km
        >>> x.value_in(units.m)
        1000.0
        
        """
        return self.conversion_factor_from(unit)
        
    def __repr__(self):
        return 'unit<'+str(self)+'>'
        
    def combine_bases(self, base1, base2):
        result = [None] * 7
        for n1, unit1 in base1:
            found = False
            for n2, unit2 in base2:
                if unit1 == unit2:
                    base2 = filter(lambda x : x[1] != unit1, base2)
                    found = True
                    result[unit1.index] = (n1, n2, unit1)
                    break
            if not found:
                result[unit1.index] =  n1, 0, unit1
        for n2, unit2 in base2:
                result[unit2.index] = 0, n2, unit2
        for x in result:
            if not x is None:
                yield x
                
    def has_same_base_as(self, other):
        """Detrmine if the base of other is the same as the
        base of self.
        
        :argument other: unit to compare base to
        :result: True, if bases are compatiple.
        
        >>> from amuse.support.units import units
        >>> mps = units.m / units.s
        >>> kph = units.km / units.hour
        >>> mps.has_same_base_as(kph)
        True
        >>> mps.has_same_base_as(units.km)
        False
        
        """
        return other.base == self.base
            
    def base_unit(self):
        if not self.base:
            return none_unit('none', 'none') # * self.factor
        
        unit = 1
        for n, base in self.base:
            if n == 1:
                unit = unit*base 
            else:
                unit = unit*(base ** n)        
        return unit
    
    def is_non_numeric(self):
        return False

    def components(self):
        return [(self,1)]

    def simplify(self):
        components=self.components()
        dic={}
        for x,i in components:
          if dic.has_key(x):
            dic[x]=dic[x]+i
          else:
            dic[x]=i
        dic=sorted(dic.items(),cmp = lambda x,y: y[1]-x[1])
        simpl=1
        for x in dic:
          if(x[1]!=0): simpl=simpl*x[0]**x[1]
        if simpl==1:
          simpl=none_unit('none','()')
        return self.conversion_factor_from(simpl)*simpl    
            
class base_unit(unit):
    """
    base_unit objects are  orthogonal, indivisable units 
    of a sytem of units.
    
    A system of units contains a set of base units
    
    :argument quantity: name of the base quantity, for example *length*
    :argument name: name of the unit, for example *meter*
    :argument symbol: symbol of the unit, for example *m*
    :argument system: system of units object
    
    >>> cgs = system("cgs")
    >>> cm = base_unit("length", "centimetre", "cm", cgs)
    >>> cm
    unit<cm>
    """
    def __init__(self, quantiy, name, symbol, system):
        self.quantity = quantiy
        self.name = name
        self.symbol = symbol
        system.add_base(self)
        
    def __str__(self):
        return self.symbol
    
    @property
    def factor(self):
        """
        The multiplication factor of a unit.
        For example, factor is 1000 for km. 
        """
        return 1
        
    @late
    def base(self):
        """
        The base represented as a list of tuples.
        Each tuple consists of an power and a unit.
        """
        return ((1,self),)
        
class none_unit(unit):
    def __init__(self, name,  symbol='()'):
        self.name = name
        self.symbol = symbol
        
    def __str__(self):
        return self.symbol
        
    def components(self):
        return []
    
    @late
    def factor(self):
        return 1
        
    @late
    def base(self):
        return ()
                    
class named_unit(unit):
    """
    A named_unit object defines an alias for another
    unit. When printing a named_unit, the symbol
    is shown and not the unit parts. For all other
    operations the named_units works exactly like
    the aliased unit.
    
    :argument name: Long name or description of the unit
    :argument symbol: Short name to show when printing units
        or quantities
    :argument unit: The unit to alias
    
    >>> from amuse.support.units import si
    >>> 60 * si.s
    unit<60 * s>
    >>> minute = named_unit("minute","min", 60*si.s)
    >>> minute
    unit<min>
    >>> (20 | (60 * si.s)).as_quantity_in(minute)
    quantity<20.0 min>
    """
    def __init__(self, name, symbol, unit):
        self.name = name
        self.symbol = symbol
        self.local_unit = unit
        
    def __str__(self):
        return self.symbol
    
    @late
    def factor(self):
        return self.local_unit.factor
        
    @late
    def base(self):
        return self.local_unit.base
        

class derived_unit(unit):
    """
    Abstract base class of derived units. New units
    can be derived from base_units. Each operation on
    a unit creates a new derived_unit.
    """
    pass
    
class factor_unit(derived_unit):
    """
    A factor_unit object defines a unit multiplied by
    a number. Do not call this method directly,
    factor_unit objects are supposed to be created by
    multiplying a number with a unit.
    
    :argument unit: The unit to derive from.
    :argument factor: The multiplication factor.
    
    >>> from amuse.support.units import si
    >>> minute = 60.0 * si.s
    >>> minute.as_quantity_in(si.s)
    quantity<60.0 s>
    >>> hour = 60.0 * minute
    >>> hour
    unit<60.0 * 60.0 * s>
    >>> hour.as_quantity_in(si.s)
    quantity<3600.0 s>
    
    """
    def __init__(self, factor , unit, name = '', symbol = ''):
        self.name = name
        self.symbol = symbol        
        self.local_factor = factor
        self.local_unit = unit
        
    def __str__(self):
        if self.symbol is '':
           return str(self.local_factor) + ' * ' + str(self.local_unit)
        return self.symbol 

    def __rmul__(self, other):
        fac=other*self.local_factor
        if fac == 1.0:
            return self.local_unit
        else:
            return factor_unit(other*self.local_factor, self.local_unit)

    def __pow__(self,other):
       return self.local_factor**other*self.local_unit**other
    
    def components(self):
       return self.local_unit.components()

    @late
    def factor(self):
        return self.local_factor * self.local_unit.factor
        
    @late
    def base(self):
        return self.local_unit.base
        
class mul_unit(derived_unit):
    """
    A mul_unit object defines a unit multiplied by
    another unit. Do not call this method directly,
    mul_unit objects are supposed to be created by
    multiplying units.
    
    :argument left_hand: Left hand side of the multiplication.
    :argument right_hand: Right hand side of the multiplication.
    
    >>> from amuse.support.units import si
    >>> area = si.m * si.m
    >>> area
    unit<m * m>
    >>> hectare = (100 * si.m) * (100 * si.m)
    >>> hectare.as_quantity_in(area)
    quantity<10000.0 m * m>
    
    """
    def __init__(self, left_hand , right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand
        
    def __str__(self):
        return str(self.left_hand) + ' * ' + str(self.right_hand) 
    
    def __pow__(self,other):
       return self.left_hand**other*self.right_hand**other

    def components(self):
       return self.left_hand.components()+self.right_hand.components()
        
    @late
    def factor(self):
        return self.left_hand.factor * self.right_hand.factor
   
    @late
    def base(self):
        return tuple(
            filter(lambda x: x[0] != 0,
            map(lambda x: (x[0] + x[1], x[2]),self.combine_bases(self.left_hand.base, self.right_hand.base))))
        
class pow_unit(derived_unit):
    """
    A pow_unit object defines a unit as
    another unit to a specified powe. 
    
    Do not call this method directly,
    pow_unit objects are supposed to be created by
    taking powers of units.
    
    :argument power: Power of the unit
    :argument unit: The unit to derive from
    
    >>> from amuse.support.units import si
    >>> area = si.m**2
    >>> area
    unit<m**2>
    >>> area.as_quantity_in(si.m * si.m)
    quantity<1.0 m * m>
    >>> hectare = (100 * si.m) ** 2
    >>> hectare.as_quantity_in(area)
    quantity<10000.0 m**2>
    
    """
    def __init__(self, power , unit):
        self.power = power
        self.local_unit = unit
        
    def __str__(self):
        if isinstance(self.local_unit, derived_unit):
          return '('+str(self.local_unit)+')' + '**' + str(self.power)
        return str(self.local_unit)+ '**' + str(self.power)
    
    def __pow__(self,other):
        power=self.power*other
        if(power==1):
          return self.local_unit
        else:
          return pow_unit(self.power*other,self.local_unit)
         
    def components(self):
       return map(lambda x: (x[0],x[1]*self.power),self.local_unit.components())

    @late
    def base(self):
        return tuple(
            filter(lambda x: x[0] != 0,
            map(lambda x : (x[0] * self.power, x[1]), self.local_unit.base)))
        
    @late
    def factor(self):
        return self.local_unit.factor ** self.power
        
class div_unit(derived_unit):
    """
    A div_unit object defines a unit multiplied by
    another unit. Do not call this method directly,
    div_unit objects are supposed to be created by
    dividing units.
    
    :argument left_hand: Left hand side of the multiplication.
    :argument right_hand: Right hand side of the multiplication.
    
    >>> from amuse.support.units import si
    >>> speed = si.m / si.s
    >>> speed
    unit<m / s>
    >>> speed_with_powers = si.m * si.s ** -1
    >>> speed.as_quantity_in(speed_with_powers)
    quantity<1.0 m * s**-1>
    """
    def __init__(self, left_hand , right_hand):
        self.left_hand = left_hand
        self.right_hand = right_hand
        
    def __pow__(self,other):
       return self.left_hand**other*self.right_hand**-other

    def __str__(self):
        if isinstance(self.right_hand, derived_unit):
          if not isinstance(self.right_hand, pow_unit):
            return str(self.left_hand) + ' / ' + '('+str(self.right_hand)+')'
        return str(self.left_hand) + ' / ' +str(self.right_hand)
        
    def components(self):
       return self.left_hand.components()+ \
         map(lambda x: (x[0],-x[1]),self.right_hand.components())

    @late
    def factor(self):
        return  self.left_hand.factor * 1.0  / self.right_hand.factor
        
    @late
    def base(self):
        return tuple(
                    filter(lambda x: x[0] != 0,
                    map(lambda x: (x[0] - x[1], x[2]),
                        self.combine_bases(self.left_hand.base, self.right_hand.base))))
        

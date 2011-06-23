from __future__ import division

class Quantity(object):
    __slots__ = ['unit']
    
    def __init__(self, unit):
        self.unit = unit
        
    def is_quantity(self):
        """
        True for all quantities.
        """
        return True
        
    def is_scalar(self):
        """
        True for scalar quantities.
        """
        return False
        
    def is_vector(self):
        """
        True for vector quantities.
        """
        return False
        
    def __repr__(self):
        return 'quantity<'+str(self)+'>'

    def __add__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.new(self.number + other_in_my_units.number , self.unit)
        
    def __sub__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.new(self.number - other_in_my_units.number , self.unit)
        
    def __ror__(self, other):
        return self.__rmul__(other)

    def __or__(self, other):
        return other|self

    def __mul__(self, other):
        if hasattr(other, "unit"):
            return self.new(self.number * other.number , (self.unit * other.unit).simplify())
#            return self.new(self.number * other.number , (self.unit * other.unit).to_simple_form())
        else:
            return self.new(self.number * other , self.unit)
            
    def __pow__(self, other):
        return self.new(self.number ** other , self.unit ** other)
        
    def __rmul__(self, other):
       return self.__mul__(other)
              
    def __div__(self, other):
        if hasattr(other, "unit"):
            return self.new(self.number / other.number , (self.unit / other.unit).simplify())
#            return self.new(self.number / other.number , (self.unit / other.unit).to_simple_form())
        else:
            return self.new(self.number / other , self.unit)
    
    def __rdiv__(self, other):
            return self.new(other / self.number , (1 / self.unit))

    __truediv__=__div__
            
    def in_(self, x):
        return self.as_quantity_in(x)
    def convert_to(self, x):
        return self.as_quantity_in(x)
            
    def in_base(self):
        unit=self.unit.base_unit()
        return self.as_quantity_in(unit)
    
    def as_quantity_in(self, another_unit): 
        value_of_unit_in_another_unit = self.unit.as_quantity_in(another_unit)
        return self.new(self.number * value_of_unit_in_another_unit.number, another_unit)

    def value_in(self, unit):
        value_of_unit_in_another_unit = self.unit.value_in(unit)
        return self.number * value_of_unit_in_another_unit

    def simplify(self):
        simpl=self.in_(self.unit.simplify())
        if not self.unit.base:
            return simpl
        if hasattr(simpl.unit,"local_factor"):
          simpl=self.new(simpl.number*simpl.unit.local_factor,simpl.unit.local_unit)
        return simpl  

    def to_simple_form(self):
        simpl=self.in_(self.unit.to_simple_form())
        if not self.unit.base:
            return simpl
        if hasattr(simpl.unit,"local_factor"):
          simpl=self.new(simpl.number*simpl.unit.local_factor,simpl.unit.local_unit)
        return simpl  


class ScalarQuantity(Quantity):
    __slots__ = ['number']
    
    def __init__(self, number, unit):
        Quantity.__init__(self, unit)
        self.number = number
                  
    def is_scalar(self):
        return True
      
    def __str__(self):
        unit_str = str(self.unit)
        if unit_str:
            return str(self.number) + ' | ' + unit_str
        else:
            return str(self.number)
                                
    def __lt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number < other_in_my_units.number
        
    def __gt__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number > other_in_my_units.number
        
    def __eq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number == other_in_my_units.number
        
    def __neq__(self, other):
        other_in_my_units = other.as_quantity_in(self.unit)
        return self.number != other_in_my_units.number

    def __le__(self, other):
        return self.__eq__(other) or self.__lt__(other)

    def __ge__(self, other):
        return self.__eq__(other) or self.__gt__(other)
            
    def copy(self):
        return self.new(self.number, self.unit)

    def to_unit(self):
        in_base=self.in_base()
        return self.number * self.unit

    def new(self,value,unit):
        if not unit.base:
            return value*unit.factor
        return ScalarQuantity(value, unit)

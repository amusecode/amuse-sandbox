from __future__ import division
import mpmath as math
from units import *

def sqrt(x):
  if hasattr(x,'unit'):
    return math.sqrt(x.number) | x.unit**0.5 
  return math.sqrt(x)  

cos=math.cos
sin=math.sin
tan=math.tan
cosh=math.cosh
sinh=math.sinh
tanh=math.tanh
acos=math.acos
asin=math.asin
atan=math.atan
acosh=math.acosh
asinh=math.asinh
atanh=math.atanh
atan2=math.atan2
cot=math.cot
coth=math.coth
j=math.j

log=math.log
log10=math.log10
fac=math.fac

e=math.e
pi=math.pi

rnd=math.rand

#unit functions:
def in_base(x):
  if hasattr(x,'unit'):
    return x.in_base()
  return x

def strip(x):
  if hasattr(x,'unit'):
    return x.number
  return x
  
def clean(x):
  if hasattr(x,'unit'):
    return x.simplify()
  return x

def simplify(x):
  if hasattr(x,'unit'):
    return x.to_simple_form()
  return x

def constants():
  return filter(lambda x: isinstance(globals()[x], PhysicalConstant),globals().keys())

def siunits():
  return filter(lambda x: isinstance(globals()[x], core.base_unit),globals().keys())

def named_units():
  return filter(lambda x: isinstance(globals()[x], core.named_unit),globals().keys())

def si_prefixes():
  return ['deca(','hecto(','kilo(','mega(','giga(','terra(','peta(', \
      'exa(','zetta(','yotta(','deci(','centi(','milli(','micro(','nano(','pico(',\
      'femto(','atto(','zepto(','yocto(']

class calcengine(object):

  def __init__(self):
    self.ans=0
    self.express=[]
    self.i_express=-1

  def pushback_result(self,result):
    self.ans=result    

  def pushback_express(self,express):
    self.express.append(express)    
    self.i_express=len(self.express)-1

  def forward_express(self):
    if(self.i_express<len(self.express)-1): self.i_express+=1
    return self.express[self.i_express]

  def backward_express(self):
    if(self.i_express>0): self.i_express-=1
    return self.express[self.i_express]

  def calculate(self,express):
    return eval(express.replace("ans","self.ans"))


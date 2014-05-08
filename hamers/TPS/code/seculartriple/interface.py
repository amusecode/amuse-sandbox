from amuse.community import *
from amuse.units import units
cm = 1e-2*units.m
g = 1e-3*units.kg

class SecularTripleInterface(CodeInterface):
    include_headers = ['src/code.h']

    def __init__(self):
         CodeInterface.__init__(self)

#    @legacy_function
#    def testf():
#        function = LegacyFunctionSpecification()
#        function.addParameter('x', 'd', function.IN)
#        function.addParameter('y', 'd', function.IN)
#        function.result_type = 'd'
#        return function

    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('R1', dtype='float64', direction=function.IN)
        function.addParameter('R2', dtype='float64', direction=function.IN)
        function.addParameter('a1', dtype='float64', direction=function.IN)
        function.addParameter('a2', dtype='float64', direction=function.IN)
        function.addParameter('e1', dtype='float64', direction=function.IN)
        function.addParameter('e2', dtype='float64', direction=function.IN)
        function.addParameter('itot', dtype='float64', direction=function.IN)
        function.addParameter('g1', dtype='float64', direction=function.IN)
        function.addParameter('g2', dtype='float64', direction=function.IN)
        function.addParameter('t', dtype='float64', direction=function.IN)        
        function.addParameter('dt', dtype='float64', direction=function.IN)        
        function.addParameter('m1_out', dtype='float64', direction=function.OUT) 
        function.addParameter('m2_out', dtype='float64', direction=function.OUT) 
        function.addParameter('m3_out', dtype='float64', direction=function.OUT) 
        function.addParameter('R1_out', dtype='float64', direction=function.OUT) 
        function.addParameter('R2_out', dtype='float64', direction=function.OUT) 
        function.addParameter('a1_out', dtype='float64', direction=function.OUT) 
        function.addParameter('a2_out', dtype='float64', direction=function.OUT) 
        function.addParameter('e1_out', dtype='float64', direction=function.OUT) 
        function.addParameter('e2_out', dtype='float64', direction=function.OUT) 
        function.addParameter('itot_out', dtype='float64', direction=function.OUT) 
        function.addParameter('g1_out', dtype='float64', direction=function.OUT) 
        function.addParameter('g2_out', dtype='float64', direction=function.OUT) 
        function.addParameter('t_out', dtype='float64', direction=function.OUT)
        function.addParameter('out_flag', dtype='int32', direction=function.OUT)
        function.addParameter('error_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function


    @legacy_function
    def get_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('relative_tolerance', dtype='float64',direction=function.OUT,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('relative_tolerance', dtype='float64',direction=function.IN,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_f_1PN_1():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_1', dtype='float64',direction=function.OUT,description = "Inner binary 1PN binary multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_oct():
        function = LegacyFunctionSpecification()
        function.addParameter('f_oct', dtype='float64',direction=function.IN,description = "Octupole term multiplication factor")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_f_oct():
        function = LegacyFunctionSpecification()
        function.addParameter('f_oct', dtype='float64',direction=function.OUT,description = "Octupole term multiplication factor")
        function.result_type = 'int32'
        return function        

    @legacy_function
    def set_f_1PN_1():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_1', dtype='float64',direction=function.IN,description = "Inner binary 1PN binary multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_1PN_2():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_2', dtype='float64',direction=function.OUT,description = "Outer binary 1PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_1PN_2():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_2', dtype='float64',direction=function.IN,description = "Outer binary 1PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_1PN_12():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_12', dtype='float64',direction=function.OUT,description = "1PN interaction multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_1PN_12():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_12', dtype='float64',direction=function.IN,description = "1PN interaction multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_25PN_1():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_1', dtype='float64',direction=function.OUT,description = "Inner binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_25PN_1():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_1', dtype='float64',direction=function.IN,description = "Inner binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_25PN_2():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_2', dtype='float64',direction=function.OUT,description = "Outer binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_25PN_2():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_2', dtype='float64',direction=function.IN,description = "Outer binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function


class SecularTriple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularTripleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_time",
            "set_time",
            "time",
            "model time",
            default_value = 0.0 | units.Myr
        )
        object.add_method_parameter(
            "get_relative_tolerance",
            "set_relative_tolerance",
            "relative_tolerance",
            "Relative tolerance, default 1e-10",
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_oct",
            "set_f_oct",
            "f_oct",
            "Octupole term multiplication factor", 
            default_value = 1
        )
                    
        object.add_method_parameter(
            "get_f_1PN_1",
            "set_f_1PN_1",
            "f_1PN_1",
            "Inner binary 1PN multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_1PN_2",
            "set_f_1PN_2",
            "f_1PN_2",
            "Outer binary 1PN multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_1PN_12",
            "set_f_1PN_12",
            "f_1PN_12",
            "1PN interaction multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_25PN_1",
            "set_f_25PN_1",
            "f_25PN_1",
            "Inner binary 2.5PN multiplication factor", 
            default_value = 0
        )
        
        object.add_method_parameter(
            "get_f_25PN_2",
            "set_f_25PN_2",
            "f_25PN_2",
            "Outer binary 2.5PN binary multiplication factor", 
            default_value = 0
        )
    
    def define_methods(self, object):

        object.add_method(
            "evolve",
            (
                g,
                g,
                g,
                cm,
                cm,
                cm,
                cm,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                units.s,
                units.s,
            ),
            (
                g,
                g,
                g,
                cm,
                cm,
                cm,
                cm,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                object.NO_UNIT,
                units.s,
                object.NO_UNIT,
                object.NO_UNIT,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_time",
            (),
            (units.Myr, object.ERROR_CODE,)
        )
        object.add_method(
            "set_time",
            (units.Myr, ),
            (object.ERROR_CODE,)
        )       
        object.add_method(
            "get_relative_tolerance",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_relative_tolerance",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_f_1PN_1",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_1PN_1",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_f_oct",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_oct",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )        

        object.add_method(
            "get_f_1PN_2",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_2",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_1PN_12",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_12",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_1",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_1",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_2",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_2",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

    def before_get_parameter(self):
        """
        Called everytime just before a parameter is retrieved in using::
            instance.parameter.name
        """
        pass
        
    def before_set_parameter(self):
        """
        Called everytime just before a parameter is updated in using::
            instance.parameter.name = newvalue
        """
        pass

    def define_particle_sets(self,object):
        object.define_inmemory_set('binaries')
        
    def evolve_model(self,end_time):
        if end_time is None:
            print 'Please specify end time!'
            return
        
        binaries = self.binaries
        m1 = binaries[0].child1.mass
        m2 = binaries[0].child2.mass
        m3 = binaries[1].child1.mass
        R1 = binaries[0].child1.radius
        R2 = binaries[0].child2.radius
        R3 = binaries[1].child1.radius
        ain = binaries[0].semimajor_axis
        ein = binaries[0].eccentricity
        aout = binaries[1].semimajor_axis
        eout = binaries[1].eccentricity
        itot = binaries[0].inclination
        gin = binaries[0].argument_of_pericenter
        gout = binaries[1].argument_of_pericenter
        
        time_step = end_time - self.model_time
#        print 'ein',ein,
        m1,m2,m3,R1,R2,ain,aout,ein,eout,itot,gin,gout,end_time_dummy,flag,error = self.evolve(m1,m2,m3,R1,R2,ain,aout,ein,eout,itot,gin,gout,self.model_time,time_step)
#        print 'ein',ein

        self.model_time = end_time
        self.binaries[0].child1.mass = m1
        self.binaries[0].child2.mass = m2
        self.binaries[1].child1.mass = m3
        self.binaries[0].child1.radius = R1
        self.binaries[0].child2.radius = R2
        self.binaries[1].child1.radius = R3
        self.binaries[0].semimajor_axis = ain
        self.binaries[0].eccentricity = ein
        self.binaries[1].semimajor_axis = aout
        self.binaries[1].eccentricity = eout
        self.binaries[0].inclination = itot
        self.binaries[0].argument_of_pericenter = gin
        self.binaries[1].argument_of_pericenter = gout
        
#        object.define_inmemory_set('binaries')
#        object.define_set('binaries', 'index_of_the_star')

#    def define_binary_sets(self,object):
#        object.define_inmemory_set('binaries')


        """
        object.define_set('particles', 'index_of_the_star')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_star')

        object.define_set('binaries', 'index_of_the_star')
        object.set_new('binaries', 'new_binary')
        object.set_delete('binaries', 'delete_binary')
        """

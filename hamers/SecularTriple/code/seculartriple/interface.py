"""
Interface for SecularTriple

Adrian Hamers 16-06-2014
"""

from amuse.community import *
from amuse.units import units
import tidal_friction_constant
cm = 1e-2*units.m
g = 1e-3*units.kg

class SecularTripleInterface(CodeInterface):
    include_headers = ['src/main_code.h']

    def __init__(self):
         CodeInterface.__init__(self)

    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('R1', dtype='float64', direction=function.IN)
        function.addParameter('R2', dtype='float64', direction=function.IN)
        function.addParameter('R3', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency1', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency2', dtype='float64', direction=function.IN)   
        function.addParameter('spin_angular_frequency3', dtype='float64', direction=function.IN)                   
        function.addParameter('AMC_star1', dtype='float64', direction=function.IN)   
        function.addParameter('AMC_star2', dtype='float64', direction=function.IN)   
        function.addParameter('AMC_star3', dtype='float64', direction=function.IN)           
        function.addParameter('gyration_radius_star1', dtype='float64', direction=function.IN)   
        function.addParameter('gyration_radius_star2', dtype='float64', direction=function.IN)   
        function.addParameter('gyration_radius_star3', dtype='float64', direction=function.IN)   
        function.addParameter('k_div_T_tides_star1', dtype='float64', direction=function.IN)   
        function.addParameter('k_div_T_tides_star2', dtype='float64', direction=function.IN)   
        function.addParameter('k_div_T_tides_star3', dtype='float64', direction=function.IN)   
        function.addParameter('a_in', dtype='float64', direction=function.IN)
        function.addParameter('a_out', dtype='float64', direction=function.IN)
        function.addParameter('e_in', dtype='float64', direction=function.IN)
        function.addParameter('e_out', dtype='float64', direction=function.IN)
        function.addParameter('INCL_in', dtype='float64', direction=function.IN)
        function.addParameter('INCL_out', dtype='float64', direction=function.IN)
#        function.addParameter('INCL12', dtype='float64', direction=function.IN)
        function.addParameter('AP_in', dtype='float64', direction=function.IN)
        function.addParameter('AP_out', dtype='float64', direction=function.IN)
        function.addParameter('LAN_in', dtype='float64', direction=function.IN)
        function.addParameter('LAN_out', dtype='float64', direction=function.IN)
        function.addParameter('m1dot', dtype='float64', direction=function.IN)
        function.addParameter('gamma_in', dtype='float64', direction=function.IN)
        function.addParameter('gamma_out', dtype='float64', direction=function.IN)
        function.addParameter('mu_in', dtype='float64', direction=function.IN)
        function.addParameter('mu_out', dtype='float64', direction=function.IN)
        function.addParameter('t', dtype='float64', direction=function.IN)        
        function.addParameter('dt', dtype='float64', direction=function.IN)        
        function.addParameter('m1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m3_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R3_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency1_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency2_output', dtype='float64', direction=function.OUT)
        function.addParameter('spin_angular_frequency3_output', dtype='float64', direction=function.OUT)
        function.addParameter('a_in_output', dtype='float64', direction=function.OUT) 
        function.addParameter('a_out_output', dtype='float64', direction=function.OUT) 
        function.addParameter('e_in_output', dtype='float64', direction=function.OUT) 
        function.addParameter('e_out_output', dtype='float64', direction=function.OUT) 
        function.addParameter('INCL_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('INCL_in_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('AP_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('AP_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('LAN_in_output', dtype='float64', direction=function.OUT)
        function.addParameter('LAN_out_output', dtype='float64', direction=function.OUT)
        function.addParameter('t_output', dtype='float64', direction=function.OUT)
        function.addParameter('output_flag', dtype='int32', direction=function.OUT)
        function.addParameter('error_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('equations_of_motion_specification', dtype='int32',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('equations_of_motion_specification', dtype='int32',direction=function.IN,description = "...")
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
    def set_f_quad():
        function = LegacyFunctionSpecification()
        function.addParameter('f_quad', dtype='float64',direction=function.IN,description = "Quadrupole term multiplication factor")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_f_quad():
        function = LegacyFunctionSpecification()
        function.addParameter('f_quad', dtype='float64',direction=function.OUT,description = "Quadrupole term multiplication factor")
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
    def set_f_mass_transfer():
        function = LegacyFunctionSpecification()
        function.addParameter('f_mass_transfer', dtype='float64',direction=function.IN,description = "Mass transfer multiplication factor")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_f_mass_transfer():
        function = LegacyFunctionSpecification()
        function.addParameter('f_mass_transfer', dtype='float64',direction=function.OUT,description = "Mass transfer multiplication factor")
        function.result_type = 'int32'
        return function        

    @legacy_function
    def set_f_tides():
        function = LegacyFunctionSpecification()
        function.addParameter('f_tides', dtype='float64',direction=function.IN,description = "Tides multiplication factor")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_f_tides():
        function = LegacyFunctionSpecification()
        function.addParameter('f_tides', dtype='float64',direction=function.OUT,description = "Tides multiplication factor")
        function.result_type = 'int32'
        return function   

    @legacy_function
    def get_f_1PN_in():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_in', dtype='float64',direction=function.OUT,description = "Outer binary 1PN multiplication factor")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_f_1PN_in():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_in', dtype='float64',direction=function.IN,description = "Inner binary 1PN binary multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_1PN_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_out', dtype='float64',direction=function.OUT,description = "Outer binary 1PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_1PN_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_out', dtype='float64',direction=function.IN,description = "Outer binary 1PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_1PN_in_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_in_out', dtype='float64',direction=function.OUT,description = "1PN interaction multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_1PN_in_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_1PN_in_out', dtype='float64',direction=function.IN,description = "1PN interaction multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_25PN_in():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_in', dtype='float64',direction=function.OUT,description = "Inner binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_25PN_in():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_in', dtype='float64',direction=function.IN,description = "Inner binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_f_25PN_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_out', dtype='float64',direction=function.OUT,description = "Outer binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_f_25PN_out():
        function = LegacyFunctionSpecification()
        function.addParameter('f_25PN_out', dtype='float64',direction=function.IN,description = "Outer binary 2.5PN multiplication factor")
        function.result_type = 'int32'
        return function


class SecularTriple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularTripleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr

    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_equations_of_motion_specification",
            "set_equations_of_motion_specification",
            "equations_of_motion_specification",
            "..",
            default_value = 0
        )
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
            "get_f_quad",
            "set_f_quad",
            "f_quad",
            "Quadrupole term multiplication factor", 
            default_value = 1
        )
        object.add_method_parameter(
            "get_f_oct",
            "set_f_oct",
            "f_oct",
            "Octupole term multiplication factor", 
            default_value = 1
        )
        object.add_method_parameter(
            "get_f_mass_transfer",
            "set_f_mass_transfer",
            "f_mass_transfer",
            "Mass transfer multiplication factor", 
            default_value = 1
        )   
        object.add_method_parameter(
            "get_f_tides",
            "set_f_tides",
            "f_tides",
            "Tides multiplication factor", 
            default_value = 1
        )                           
        object.add_method_parameter(
            "get_f_1PN_in",
            "set_f_1PN_in",
            "f_1PN_in",
            "Inner binary 1PN multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_1PN_out",
            "set_f_1PN_out",
            "f_1PN_out",
            "Outer binary 1PN multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_1PN_in_out",
            "set_f_1PN_in_out",
            "f_1PN_in_out",
            "1PN interaction multiplication factor", 
            default_value = 0
        )

        object.add_method_parameter(
            "get_f_25PN_in",
            "set_f_25PN_in",
            "f_25PN_in",
            "Inner binary 2.5PN multiplication factor", 
            default_value = 0
        )
        
        object.add_method_parameter(
            "get_f_25PN_out",
            "set_f_25PN_out",
            "f_25PN_out",
            "Outer binary 2.5PN binary multiplication factor", 
            default_value = 0
        )
    
    def define_methods(self, object):

        object.add_method(
            "evolve",
            (
                g,                          ### m1
                g,                          ### m2
                g,                          ### m3
                cm,                         ### R1
                cm,                         ### R2
                cm,                         ### R3
                1.0/units.s,                ### spin_angular_frequency1
                1.0/units.s,                ### spin_angular_frequency2
                1.0/units.s,                ### spin_angular_frequency3                                
                object.NO_UNIT,             ### AMC_star1
                object.NO_UNIT,             ### AMC_star2
                object.NO_UNIT,             ### AMC_star3
                object.NO_UNIT,             ### gyration_radius_star1
                object.NO_UNIT,             ### gyration_radius_star2
                object.NO_UNIT,             ### gyration_radius_star3
                1.0/units.s,                ### k_div_T_tides_star1
                1.0/units.s,                ### k_div_T_tides_star2
                1.0/units.s,                ### k_div_T_tides_star3                                
                cm,                         ### a_in
                cm,                         ### a_out
                object.NO_UNIT,             ### e_in
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### INCL_in
                object.NO_UNIT,             ### INCL_out                
#                object.NO_UNIT,             ### i12                
                object.NO_UNIT,             ### AP_in
                object.NO_UNIT,             ### AP_out
                object.NO_UNIT,             ### LAN_in
                object.NO_UNIT,             ### LAN_out
                g/units.s,                  ### m1dot       
                object.NO_UNIT,             ### gamma_in
                object.NO_UNIT,             ### gamma_out
                object.NO_UNIT,             ### mu_in
                object.NO_UNIT,             ### mu_out
                units.s,                    ### t
                units.s,                    ### dt
            ),
            (
                g,                          ### m1_output
                g,                          ### m2_output
                g,                          ### m3_output
                cm,                         ### R1_output
                cm,                         ### R2_output
                cm,                         ### R3_output
                1.0/units.s,                ### spin_angular_frequency1_output
                1.0/units.s,                ### spin_angular_frequency2_output
                1.0/units.s,                ### spin_angular_frequency3_output                         
                cm,                         ### a_in_output
                cm,                         ### a_out_output
                object.NO_UNIT,             ### e_in_output
                object.NO_UNIT,             ### e_out_output
                object.NO_UNIT,             ### INCL_in_output
                object.NO_UNIT,             ### INCL_out_output
                object.NO_UNIT,             ### INCL_in_out_output
                object.NO_UNIT,             ### AP_in_output
                object.NO_UNIT,             ### AP_out_output                                
                object.NO_UNIT,             ### LAN_in_output
                object.NO_UNIT,             ### LAN_out_output                                
                units.s,                    ### t_output
                object.NO_UNIT,             ### output_flag
                object.NO_UNIT,           ### error_flag
                object.ERROR_CODE,
            )
        )
        
        object.add_method(
            "get_equations_of_motion_specification",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_equations_of_motion_specification",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
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
            "get_f_quad",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_quad",
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
            "get_f_mass_tranfer",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_mass_tranfer",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )   
        object.add_method(
            "get_f_tides",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_tides",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_f_1PN_in",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_f_1PN_in",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )        
        object.add_method(
            "get_f_1PN_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_out",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_1PN_in_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_1PN_in_out",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_in",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_in",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )

        object.add_method(
            "get_f_25PN_out",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        
        object.add_method(
            "set_f_25PN_out",
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

        ### extract data from binaries ###
        binaries = self.binaries
        star1 = binaries[0].child1
        star2 = binaries[0].child2
        star3 = binaries[1].child1

        m1 = star1.mass
        m2 = star2.mass
        m3 = star3.mass

        a_in = binaries[0].semimajor_axis
        e_in = binaries[0].eccentricity
        a_out = binaries[1].semimajor_axis
        e_out = binaries[1].eccentricity
        
        INCL_in = binaries[0].inclination
        INCL_out = binaries[1].inclination
        AP_in = binaries[0].argument_of_pericenter
        AP_out = binaries[1].argument_of_pericenter
        LAN_in = binaries[0].longitude_of_ascending_node
        LAN_out = binaries[1].longitude_of_ascending_node        
        
        do_tides = False
        do_mass_transfer = False
        
        if self.parameters.f_mass_transfer != 0.0: do_mass_transfer = True
        if self.parameters.f_tides!=0.0: do_tides = True        

        ### following quantities are not used if mass transfer/tides are not taken into account, but need to pass them to evolve() ###
        R1 = R2 = R3 = 0.0 | units.RSun
        m1dot = 0.0 | units.MSun/units.yr
        gamma_in = gamma_out = mu_in = mu_out=0.0
        spin_angular_frequency1=spin_angular_frequency2=spin_angular_frequency3 = 0.0 | 1.0/units.s
        AMC_star1 = AMC_star2 = AMC_star3 = 0.0
        gyration_radius_star1 = gyration_radius_star2 = gyration_radius_star3 = 0.0
        k_div_T_tides_star1 = k_div_T_tides_star2 = k_div_T_tides_star3 = 0.0 | 1.0/units.s

        if do_mass_transfer == True:
            try:
                m1dot = star1.mass_transfer_rate
                gamma_in = binaries[0].mass_transfer_accretion_parameter
                gamma_out = binaries[1].mass_transfer_accretion_parameter
                mu_in = binaries[0].mass_transfer_angular_momentum_loss_parameter
                mu_out = binaries[1].mass_transfer_angular_momentum_loss_parameter                
            except AttributeError:
                print "More attributes required for mass transfer"
                return

        if do_tides == True:
            try:
                stellar_type1 = star1.stellar_type
                stellar_type2 = star2.stellar_type
                stellar_type3 = star3.stellar_type                
                m1_envelope = star1.envelope_mass
                m2_envelope = star2.envelope_mass
                m3_envelope = star3.envelope_mass
                R1 = star1.radius
                R2 = star2.radius
                R3 = star3.radius
                R1_envelope = star1.envelope_radius
                R2_envelope = star2.envelope_radius
                R3_envelope = star3.envelope_radius                
                luminosity_star1 = star1.luminosity
                luminosity_star2 = star2.luminosity
                luminosity_star3 = star3.luminosity
                AMC_star1 = star1.apsidal_motion_constant
                AMC_star2 = star2.apsidal_motion_constant
                AMC_star3 = star3.apsidal_motion_constant
                gyration_radius_star1 = star1.gyration_radius
                gyration_radius_star2 = star2.gyration_radius
                gyration_radius_star3 = star3.gyration_radius
                spin_angular_frequency1 = star1.spin_angular_frequency
                spin_angular_frequency2 = star2.spin_angular_frequency
                spin_angular_frequency3 = star3.spin_angular_frequency
                
                k_div_T_tides_star1 = tidal_friction_constant.tidal_friction_constant(stellar_type1,m1,m2,a_in,R1,m1_envelope,R1_envelope,luminosity_star1,spin_angular_frequency1,gyration_radius_star1)
                k_div_T_tides_star2 = tidal_friction_constant.tidal_friction_constant(stellar_type1,m2,m1,a_in,R2,m2_envelope,R2_envelope,luminosity_star2,spin_angular_frequency2,gyration_radius_star2)
                k_div_T_tides_star3 = tidal_friction_constant.tidal_friction_constant(stellar_type1,m3,m1+m2,a_out,R3,m3_envelope,R3_envelope,luminosity_star3,spin_angular_frequency3,gyration_radius_star3)              
            except AttributeError:
                print "More attributes required for tides"
                return

        ### solve system of ODEs ###
        time_step = end_time - self.model_time
        m1,m2,m3,R1,R2,R3,spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,INCL_in_out,AP_in,AP_out,LAN_in,LAN_out,end_time_dummy,flag,error = self.evolve(
            m1,m2,m3,
            R1,R2,R3,
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3,
            AMC_star1,AMC_star2,AMC_star3,
            gyration_radius_star1,gyration_radius_star2,gyration_radius_star3,
            k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3,
            a_in,a_out,
            e_in,e_out,
            INCL_in,INCL_out,
            AP_in,AP_out,LAN_in,LAN_out,
            m1dot,
            gamma_in,gamma_out,mu_in,mu_out,
            self.model_time,time_step)

        ### update binaries ###
        self.model_time = end_time
#        self.binaries[0].child1.mass = m1
#        self.binaries[0].child2.mass = m2
#        self.binaries[1].child1.mass = m3
#        self.binaries[0].child1.radius = R1
#        self.binaries[0].child2.radius = R2
#        self.binaries[1].child1.radius = R3
        if do_tides == True:
            self.binaries[0].child1.spin_angular_frequency = spin_angular_frequency1
            self.binaries[0].child2.spin_angular_frequency = spin_angular_frequency2
            self.binaries[1].child1.spin_angular_frequency = spin_angular_frequency3

        self.binaries[0].semimajor_axis = a_in
        self.binaries[0].eccentricity = e_in
        self.binaries[1].semimajor_axis = a_out
        self.binaries[1].eccentricity = e_out
        self.binaries[0].inclination = INCL_in
        self.binaries[1].inclination = INCL_out
        self.binaries[0].mutual_inclination = INCL_in_out
        self.binaries[1].mutual_inclination = INCL_in_out
        self.binaries[0].argument_of_pericenter = AP_in
        self.binaries[1].argument_of_pericenter = AP_out
        self.binaries[0].longitude_of_ascending_node = LAN_in
        self.binaries[1].longitude_of_ascending_node = LAN_out

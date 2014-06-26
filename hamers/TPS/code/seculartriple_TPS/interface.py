"""
Interface for SecularTriple

Adrian Hamers 24-06-2014
"""

from amuse.community import *
from amuse.units import units,constants
import tidal_friction_constant
cm = 1e-2*units.m
g = 1e-3*units.kg

class SecularTripleInterface(CodeInterface):
    include_headers = ['src/main_code.h','src/ODE_system.h']

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
        function.addParameter('AP_in', dtype='float64', direction=function.IN)
        function.addParameter('AP_out', dtype='float64', direction=function.IN)
        function.addParameter('LAN_in', dtype='float64', direction=function.IN)
        function.addParameter('LAN_out', dtype='float64', direction=function.IN)
        function.addParameter('star1_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('star2_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('star3_is_donor', dtype='bool', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star1', dtype='float64', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star2', dtype='float64', direction=function.IN)
        function.addParameter('wind_mass_loss_rate_star3', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star1', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star2', dtype='float64', direction=function.IN)
        function.addParameter('time_derivative_of_radius_star3', dtype='float64', direction=function.IN)        
        function.addParameter('inner_mass_transfer_rate', dtype='float64', direction=function.IN)
        function.addParameter('outer_mass_transfer_rate', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_wind_child1_to_child2', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_wind_child2_to_child1', dtype='float64', direction=function.IN)
        function.addParameter('outer_accretion_efficiency_wind_child1_to_child2', dtype='float64', direction=function.IN)
        function.addParameter('outer_accretion_efficiency_wind_child2_to_child1', dtype='float64', direction=function.IN)
        function.addParameter('inner_accretion_efficiency_mass_transfer', dtype='float64', direction=function.IN)        
        function.addParameter('outer_accretion_efficiency_mass_transfer', dtype='float64', direction=function.IN)                      
        function.addParameter('inner_specific_AM_loss_mass_transfer', dtype='float64', direction=function.IN)                        
        function.addParameter('outer_specific_AM_loss_mass_transfer', dtype='float64', direction=function.IN)
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
    def roche_radius_pericenter_sepinsky():
        function = LegacyFunctionSpecification()
        function.addParameter('rp', dtype='float64', direction=function.IN)
        function.addParameter('q', dtype='float64', direction=function.IN)
        function.addParameter('e', dtype='float64', direction=function.IN)
        function.addParameter('f', dtype='float64', direction=function.IN)
        function.result_type = 'float64'
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
    def get_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_dynamical_stability', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_dynamical_stability', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_inner_collision', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_inner_collision', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_outer_collision', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_outer_collision', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_inner_RLOF', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_inner_RLOF', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_outer_RLOF', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('check_for_outer_RLOF', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_quadrupole_terms', dtype='bool',direction=function.IN,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_quadrupole_terms', dtype='bool',direction=function.OUT,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function         

    @legacy_function
    def set_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_octupole_terms', dtype='bool',direction=function.IN,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_octupole_terms', dtype='bool',direction=function.OUT,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function   

    @legacy_function
    def get_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_inner_terms', dtype='bool',direction=function.OUT,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_inner_terms', dtype='bool',direction=function.IN,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_outer_terms', dtype='bool',direction=function.OUT,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_outer_terms', dtype='bool',direction=function.IN,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_inner_outer_terms', dtype='bool',direction=function.OUT,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_1PN_inner_outer_terms', dtype='bool',direction=function.IN,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_25PN_inner_terms', dtype='bool',direction=function.OUT,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_25PN_inner_terms', dtype='bool',direction=function.IN,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_25PN_outer_terms', dtype='bool',direction=function.OUT,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_25PN_outer_terms', dtype='bool',direction=function.IN,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_tidal_terms', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_tidal_terms', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_tidal_terms', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_tidal_terms', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_inner_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_wind_terms', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_wind_terms', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_wind_terms', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_wind_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_wind_terms', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_wind_spin_coupling_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_wind_spin_coupling_terms', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_wind_spin_coupling_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_wind_spin_coupling_terms', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_inner_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_RLOF_terms', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_inner_RLOF_terms', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_RLOF_terms', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_RLOF_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('include_outer_RLOF_terms', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function
        
class SecularTriple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularTripleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr

    def define_parameters(self, object):
        
        object.add_method_parameter(
            "get_model_time",
            "set_model_time",
            "model_time",
            "model_time",
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
            "get_equations_of_motion_specification",
            "set_equations_of_motion_specification",
            "equations_of_motion_specification",
            "equations_of_motion_specification",
            default_value = 0
        )
        object.add_method_parameter(
            "get_check_for_dynamical_stability",
            "set_check_for_dynamical_stability",
            "check_for_dynamical_stability",
            "check_for_dynamical_stability",
            default_value = True
        )       
        object.add_method_parameter(
            "get_check_for_inner_collision",
            "set_check_for_inner_collision",
            "check_for_inner_collision",
            "check_for_inner_collision",
            default_value = True
        )          
        object.add_method_parameter(
            "get_check_for_outer_collision",
            "set_check_for_outer_collision",
            "check_for_outer_collision",
            "check_for_outer_collision",
            default_value = False
        )  
        object.add_method_parameter(
            "get_check_for_inner_RLOF",
            "set_check_for_inner_RLOF",
            "check_for_inner_RLOF",
            "check_for_inner_RLOF",
            default_value = False
        )
        object.add_method_parameter(
            "get_check_for_outer_RLOF",
            "set_check_for_outer_RLOF",
            "check_for_outer_RLOF",
            "check_for_outer_RLOF",
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_quadrupole_terms",
            "set_include_quadrupole_terms",
            "include_quadrupole_terms",
            "Whether or not to include quadrupole terms in the secular equations of motion", 
            default_value = True
        )        
        object.add_method_parameter(
            "get_include_octupole_terms",
            "set_include_octupole_terms",
            "include_octupole_terms",
            "Whether or not to include octupole terms in the secular equations of motion", 
            default_value = True
        )
        object.add_method_parameter(
            "get_include_1PN_inner_terms",
            "set_include_1PN_inner_terms",
            "include_1PN_inner_terms",
            "Include 1PN inner binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_1PN_inner_terms",
            "set_include_1PN_inner_terms",
            "include_1PN_inner_terms",
            "Include 1PN inner binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_1PN_outer_terms",
            "set_include_1PN_outer_terms",
            "include_1PN_outer_terms",
            "Include 1PN outer binary terms in the equations of motion", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_1PN_inner_outer_terms",
            "set_include_1PN_inner_outer_terms",
            "include_1PN_inner_outer_terms",
            "Include 1PN 'interaction terms in the equations of motion", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_25PN_inner_terms",
            "set_include_25PN_inner_terms",
            "include_25PN_inner_terms",
            "Include 25PN inner binary terms in the equations of motion", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_25PN_outer_terms",
            "set_include_25PN_outer_terms",
            "include_25PN_outer_terms",
            "Include 25PN outer binary terms in the equations of motion", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_inner_tidal_terms",
            "set_include_inner_tidal_terms",
            "include_inner_tidal_terms",
            "Whether or not to include the effects of tides in the inner binary system", 
            default_value = False
        )                           
        object.add_method_parameter(
            "get_include_outer_tidal_terms",
            "set_include_outer_tidal_terms",
            "include_outer_tidal_terms",
            "Whether or not to include the effects of tides in the outer binary system", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_inner_wind_terms",
            "set_include_inner_wind_terms",
            "include_inner_wind_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_outer_wind_terms",
            "set_include_outer_wind_terms",
            "include_outer_wind_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_wind_spin_coupling_terms",
            "set_include_wind_spin_coupling_terms",
            "include_wind_spin_coupling_terms",
            "..", 
            default_value = False
        )        
        object.add_method_parameter(
            "get_include_inner_RLOF_terms",
            "set_include_inner_RLOF_terms",
            "include_inner_RLOF_terms",
            "..", 
            default_value = False
        )
        object.add_method_parameter(
            "get_include_outer_RLOF_terms",
            "set_include_outer_RLOF_terms",
            "include_outer_RLOF_terms",
            "..", 
            default_value = False
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
                object.NO_UNIT,             ### AP_in
                object.NO_UNIT,             ### AP_out
                object.NO_UNIT,             ### LAN_in
                object.NO_UNIT,             ### LAN_out
                object.NO_UNIT,             ### star1_is_donor
                object.NO_UNIT,             ### star2_is_donor
                object.NO_UNIT,             ### star3_is_donor
                g/units.s,                  ### wind_mass_loss_rate_star1
                g/units.s,                  ### wind_mass_loss_rate_star2
                g/units.s,                  ### wind_mass_loss_rate_star3                                
                cm/units.s,                 ### time_derivative_of_radius_star1
                cm/units.s,                 ### time_derivative_of_radius_star2
                cm/units.s,                 ### time_derivative_of_radius_star3                                
                g/units.s,                  ### inner_mass_transfer_rate                                
                g/units.s,                  ### outer_mass_transfer_rate                                                
                object.NO_UNIT,             ### inner_accretion_efficiency_wind_child1_to_child2
                object.NO_UNIT,             ### inner_accretion_efficiency_wind_child2_to_child1
                object.NO_UNIT,             ### outer_accretion_efficiency_wind_child1_to_child2
                object.NO_UNIT,             ### outer_accretion_efficiency_wind_child2_to_child1
                object.NO_UNIT,             ### inner_accretion_efficiency_mass_transfer
                object.NO_UNIT,             ### outer_accretion_efficiency_mass_transfer
                object.NO_UNIT,             ### inner_specific_AM_loss_mass_transfer
                object.NO_UNIT,             ### outer_specific_AM_loss_mass_transfer
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
                object.NO_UNIT,             ### error_flag
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "roche_radius_pericenter_sepinsky",
            (
                cm,                         ### rp
                object.NO_UNIT,             ### q
                object.NO_UNIT,             ### e
                object.NO_UNIT,             ### f
            ),
            (
                cm,                         ### Roche radius
            )
        )

        """
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
            "get_include_quadrupole_terms",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_include_quadrupole_terms",
            (object.NO_UNIT, ),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_include_octupole_terms",
            (),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "set_include_octupole_terms",
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
        """
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
        object.define_inmemory_set('triples')

    def evolve_model(self,end_time):
        if end_time is None:
            print 'Please specify end time!'
            return

        parameters = self.parameters
        triples = self.triples
        
        for index_triple, triple in enumerate(triples):

            ### extract data from triple ###
            self.time_step = end_time - self.model_time            
            args = extract_data_from_triple(self,triple)

            ### solve system of ODEs ###
            m1,m2,m3,R1,R2,R3, \
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3, \
            a_in,a_out,e_in,e_out, \
            INCL_in,INCL_out,INCL_in_out, \
            AP_in,AP_out,LAN_in,LAN_out, \
            end_time_dummy,flag,error = self.evolve(*args)
            """
            m1,m2,m3,R1,R2,R3, \
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3, \
            a_in,a_out,e_in,e_out, \
            INCL_in,INCL_out,INCL_in_out, \
            AP_in,AP_out,LAN_in,LAN_out, \
            end_time_dummy,flag,error = self.evolve(
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
                star1_is_donor,star2_is_donor,star3_is_donor,
                wind_mass_loss_rate_star1,wind_mass_loss_rate_star2,wind_mass_loss_rate_star3,
                time_derivative_of_radius_star1,time_derivative_of_radius_star2,time_derivative_of_radius_star3,
                inner_mass_transfer_rate,outer_mass_transfer_rate,
                inner_accretion_efficiency_wind_child1_to_child2,inner_accretion_efficiency_wind_child2_to_child1,
                outer_accretion_efficiency_wind_child1_to_child2,outer_accretion_efficiency_wind_child2_to_child1,
                inner_accretion_efficiency_mass_transfer,outer_accretion_efficiency_mass_transfer,
                inner_specific_AM_loss_mass_transfer,outer_specific_AM_loss_mass_transfer,
                self.model_time,time_step)
            """
            
            ### dynamical instability ###
            triple.dynamical_instability = False
            if flag==1: 
                triple.dynamical_instability = True

            ### inner collision ###
            triple.inner_collision = False
            if flag==2: 
                triple.inner_collision = True

            ### outer collision ###
            triple.outer_collision = False
            if flag==3: 
                triple.outer_collision = True

            ### RLOF star1 ###
            triple.inner_binary.child1.is_donor = False
            if flag==4: 
                triple.inner_binary.child1.is_donor = True

            ### RLOF star2 ###
            triple.inner_binary.child2.is_donor = False
            if flag==5: 
                triple.inner_binary.child2.is_donor = True

            ### RLOF star3 ###
            triple.outer_binary.child1.is_donor = False
            if flag==6: 
                triple.outer_binary.child1.is_donor = True

            ### update model time ###
            self.model_time = end_time

            ### update triple particle ###
            if parameters.include_inner_tidal_terms == True:
                triple.inner_binary.child1.spin_angular_frequency = spin_angular_frequency1
                triple.inner_binary.child2.spin_angular_frequency = spin_angular_frequency2
            if parameters.include_outer_tidal_terms == True:
                triple.outer_binary.child1.spin_angular_frequency = spin_angular_frequency3

            if parameters.include_wind_spin_coupling_terms == True:
                triple.inner_binary.child1.spin_angular_frequency = spin_angular_frequency1
                triple.inner_binary.child2.spin_angular_frequency = spin_angular_frequency2
                triple.outer_binary.child1.spin_angular_frequency = spin_angular_frequency3
                                                
            triple.inner_binary.semimajor_axis = a_in
            triple.inner_binary.eccentricity = e_in
            triple.outer_binary.semimajor_axis = a_out
            triple.outer_binary.eccentricity = e_out
            triple.mutual_inclination = INCL_in
            triple.inner_binary.argument_of_pericenter = AP_in
            triple.outer_binary.argument_of_pericenter = AP_out
            triple.inner_binary.longitude_of_ascending_node = LAN_in
            triple.outer_binary.longitude_of_ascending_node = LAN_out

    def give_roche_radii(self,triple):
        if triple is None:
            print 'Please give triple particle'
            return

        inner_binary = triple.inner_binary
        outer_binary = triple.outer_binary

        star1 = inner_binary.child1
        star2 = inner_binary.child2
        star3 = outer_binary.child1

        m1 = star1.mass
        m2 = star2.mass
        m3 = star3.mass

        R1 = star1.radius
        R2 = star2.radius
        R3 = star3.radius
                    
        a_in = inner_binary.semimajor_axis
        e_in = inner_binary.eccentricity
        a_out = outer_binary.semimajor_axis
        e_out = outer_binary.eccentricity

        spin_angular_frequency1 = star1.spin_angular_frequency
        spin_angular_frequency2 = star2.spin_angular_frequency
        spin_angular_frequency3 = star3.spin_angular_frequency

        rp_in = a_in*(1.0-e_in)
        rp_out = a_out*(1.0-e_out)

        spin_angular_frequency_inner_orbit_periapse = numpy.sqrt( constants.G*(m1+m2)*(1.0+e_in)/(rp_in**3) ) 
        spin_angular_frequency_outer_orbit_periapse = numpy.sqrt( constants.G*(m1+m2+m3)*(1.0+e_out)/(rp_out**3) ) 

        f1 = spin_angular_frequency1/spin_angular_frequency_inner_orbit_periapse
        f2 = spin_angular_frequency2/spin_angular_frequency_inner_orbit_periapse
        f3 = spin_angular_frequency3/spin_angular_frequency_outer_orbit_periapse      
        
        R_L_star1 = self.roche_radius_pericenter_sepinsky(rp_in,m1/m2,e_in,f1)
        R_L_star2 = self.roche_radius_pericenter_sepinsky(rp_in,m2/m1,e_in,f2)
        R_L_star3 = self.roche_radius_pericenter_sepinsky(rp_out,m3/(m1+m2),e_out,f3)
                        
        if 1==0: ### test with circular & synchronous orbits: compare to Eggleton
            f1=f2=f3=1.0
            e_in=e_out=0.0
            rp_in = a_in*(1.0-e_in)
            rp_out = a_out*(1.0-e_out)
            
            print 'R_Ls',R_L_star1,R_L_star2,R_L_star3
            print 'Egg',R_L_eggleton(a_in,m1/m2),R_L_eggleton(a_in,m2/m1),R_L_eggleton(a_out,m3/(m1+m2))
            print 'ratios',R_L_star1/R_L_eggleton(a_in,m1/m2),R_L_star2/R_L_eggleton(a_in,m2/m1),R_L_star3/R_L_eggleton(a_out,m3/(m1+m2))
        
        return R_L_star1,R_L_star2,R_L_star3

def R_L_eggleton(a,q):
    q_pow_one_third = pow(q,1.0/3.0)
    q_pow_two_third = q_pow_one_third*q_pow_one_third
    return a*0.49*q_pow_two_third/(0.6*q_pow_two_third + numpy.log(1.0 + q_pow_one_third))

def extract_data_from_triple(self,triple):
    parameters = self.parameters
    
    ### general parameters ###
    inner_binary = triple.inner_binary
    outer_binary = triple.outer_binary

    star1 = inner_binary.child1
    star2 = inner_binary.child2
    star3 = outer_binary.child1

    m1 = star1.mass
    m2 = star2.mass
    m3 = star3.mass

    R1 = star1.radius
    R2 = star2.radius
    R3 = star3.radius
                    
    a_in = inner_binary.semimajor_axis
    e_in = inner_binary.eccentricity
    a_out = outer_binary.semimajor_axis
    e_out = outer_binary.eccentricity
        
    INCL_in = triple.mutual_inclination
    INCL_out = 0.0
    AP_in = inner_binary.argument_of_pericenter
    AP_out = outer_binary.argument_of_pericenter
    LAN_in = inner_binary.longitude_of_ascending_node
    LAN_out = outer_binary.longitude_of_ascending_node        

    ### mass variation parameters ###
    star1_is_donor = star2_is_donor = star3_is_donor = False
    wind_mass_loss_rate_star1 = wind_mass_loss_rate_star2 = wind_mass_loss_rate_star3 = 0.0 | units.MSun/units.yr
    inner_mass_transfer_rate = outer_mass_transfer_rate = 0.0 | units.MSun/units.yr
    inner_accretion_efficiency_wind_child1_to_child2 = inner_accretion_efficiency_wind_child2_to_child1 = 0.0
    outer_accretion_efficiency_wind_child1_to_child2 = outer_accretion_efficiency_wind_child2_to_child1 = 0.0
    inner_accretion_efficiency_mass_transfer = outer_accretion_efficiency_mass_transfer = 0.0
    inner_specific_AM_loss_mass_transfer = outer_specific_AM_loss_mass_transfer = 0.0

    if parameters.include_inner_wind_terms == True:
        try:
            wind_mass_loss_rate_star1 = star1.wind_mass_loss_rate
            wind_mass_loss_rate_star2 = star2.wind_mass_loss_rate
            inner_accretion_efficiency_wind_child1_to_child2 = inner_binary.accretion_efficiency_wind_child1_to_child2
            inner_accretion_efficiency_wind_child2_to_child1 = inner_binary.accretion_efficiency_wind_child2_to_child1
        except AttributeError:
            print "More attributes required for inner wind"
            return
    if parameters.include_outer_wind_terms == True:
        try:
            wind_mass_loss_rate_star3 = star3.wind_mass_loss_rate                    
            outer_accretion_efficiency_wind_child1_to_child2 = outer_binary.accretion_efficiency_wind_child1_to_child2
            outer_accretion_efficiency_wind_child2_to_child1 = outer_binary.accretion_efficiency_wind_child2_to_child1
        except AttributeError:
            print "More attributes required for outer wind"
            return
    if parameters.include_inner_RLOF_terms == True:
        try:
            star1_is_donor = star1.is_donor
            star2_is_donor = star2.is_donor
            inner_mass_transfer_rate = inner_binary.mass_transfer_rate
            inner_accretion_efficiency_mass_transfer = inner_binary.accretion_efficiency_mass_transfer
            inner_specific_AM_loss_mass_transfer = inner_binary.specific_AM_loss_mass_transfer
        except AttributeError:
            print "More attributes required for inner RLOF"
            return
    if parameters.include_outer_RLOF_terms == True:
        try:
            star3_is_donor = star3.is_donor            
            outer_mass_transfer_rate = outer_binary.mass_transfer_rate
            outer_accretion_efficiency_mass_transfer = outer_binary.accretion_efficiency_mass_transfer
            outer_specific_AM_loss_mass_transfer = outer_binary.specific_AM_loss_mass_transfer
        except AttributeError:
            print "More attributes required for outer RLOF"
            return

    ### RLOF checks ###
    spin_angular_frequency1 = spin_angular_frequency2 = spin_angular_frequency3 = 0.0 | 1.0/units.s
    gyration_radius_star1 = gyration_radius_star2 = gyration_radius_star3 = 0.0
    if parameters.check_for_inner_RLOF == True:
        try:
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
        except AttributeError:
            print "More attributes required for inner RLOF check"
            return
    if parameters.check_for_outer_RLOF == True:
        try:
            spin_angular_frequency3 = star3.spin_angular_frequency
        except AttributeError:
            print "More attributes required for outer RLOF check"
            return

    ### wind-spin coupling ###
    time_derivative_of_radius_star1=time_derivative_of_radius_star2=time_derivative_of_radius_star3 = 0.0 | units.RSun/units.s
    if parameters.include_wind_spin_coupling_terms == True:
        try:
            spin_angular_frequency1 = star1.spin_angular_frequency
            spin_angular_frequency2 = star2.spin_angular_frequency
            spin_angular_frequency3 = star3.spin_angular_frequency
            time_derivative_of_radius_star1 = star1.time_derivative_of_radius
            time_derivative_of_radius_star2 = star2.time_derivative_of_radius
            time_derivative_of_radius_star3 = star3.time_derivative_of_radius
            gyration_radius_star1 = star1.gyration_radius
            gyration_radius_star2 = star2.gyration_radius                    
            gyration_radius_star3 = star3.gyration_radius                    
        except AttributeError:
            print "More attributes required for wind-spin coupling terms"
            return

    ### tides ###
    AMC_star1 = AMC_star2 = AMC_star3 = 0.0
    k_div_T_tides_star1 = k_div_T_tides_star2 = k_div_T_tides_star3 = 0.0 | 1.0/units.s
    
    if parameters.include_inner_tidal_terms == True:
        try:
            stellar_type1 = star1.stellar_type
            stellar_type2 = star2.stellar_type
            stellar_type3 = star3.stellar_type                
            m1_envelope = star1.envelope_mass
            m2_envelope = star2.envelope_mass
            m3_envelope = star3.envelope_mass
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

    args = [m1,m2,m3,
        R1,R2,R3,
        spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3,
        AMC_star1,AMC_star2,AMC_star3,
        gyration_radius_star1,gyration_radius_star2,gyration_radius_star3,
        k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3,
        a_in,a_out,
        e_in,e_out,
        INCL_in,INCL_out,
        AP_in,AP_out,LAN_in,LAN_out,
        star1_is_donor,star2_is_donor,star3_is_donor,
        wind_mass_loss_rate_star1,wind_mass_loss_rate_star2,wind_mass_loss_rate_star3,
        time_derivative_of_radius_star1,time_derivative_of_radius_star2,time_derivative_of_radius_star3,
        inner_mass_transfer_rate,outer_mass_transfer_rate,
        inner_accretion_efficiency_wind_child1_to_child2,inner_accretion_efficiency_wind_child2_to_child1,
        outer_accretion_efficiency_wind_child1_to_child2,outer_accretion_efficiency_wind_child2_to_child1,
        inner_accretion_efficiency_mass_transfer,outer_accretion_efficiency_mass_transfer,
        inner_specific_AM_loss_mass_transfer,outer_specific_AM_loss_mass_transfer,
        self.model_time,self.time_step]
    return args

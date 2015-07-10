"""
Interface for SecularTriple

Adrian Hamers 2015
"""

from amuse.community import *
from amuse.units import units,constants

import tidal_friction_constant

### units used internally in the ODE solver ###
unit_l = units.AU
unit_m = units.MSun
unit_t = 1.0e6*units.yr
unit_h = unit_m*unit_l**2/unit_t ### angular momentum
unit_e = unit_m*unit_l**2/(unit_t**2) ### energy

### CVODE flags ###
CV_SUCCESS=0
CV_ROOT_RETURN=2
CV_WARNING=99
CV_TOO_MUCH_WORK=-1
CV_TOO_MUCH_ACC=-2
CV_ERR_FAILURE=-3
CV_CONV_FAILURE=-4
CV_LINIT_FAIL=-5
CV_LSETUP_FAIL=-6
CV_LSOLVE_FAIL=-7
CV_RHSFUNC_FAIL=-8
CV_FIRST_RHSFUNC_ERR=-9
CV_REPTD_RHSFUNC_ERR=-10
CV_UNREC_RHSFUNC_ERR=-11
CV_RTFUNC_FAIL=-12
CV_MEM_FAIL=-20
CV_MEM_NULL=-21
CV_ILL_INPUT=-22
CV_NO_MALLOC=-23
CV_BAD_K=-24
CV_BAD_T=-25
CV_BAD_DKY=-26
CV_TOO_CLOSE=-27


class SecularTripleInterface(CodeInterface):
    include_headers = ['src/main_code.h','src/ODE_system.h']

    def __init__(self, **options):
         CodeInterface.__init__(self, **options)
        
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
        function.addParameter('t', dtype='float64', direction=function.IN)        
        function.addParameter('global_time_step', dtype='float64', direction=function.IN)        
        function.addParameter('m1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('m3_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R1_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R2_output', dtype='float64', direction=function.OUT) 
        function.addParameter('R3_output', dtype='float64', direction=function.OUT)
        function.addParameter('tidal_E1_dot_output', dtype='float64', direction=function.OUT)   
        function.addParameter('tidal_E2_dot_output', dtype='float64', direction=function.OUT)   
        function.addParameter('tidal_E3_dot_output', dtype='float64', direction=function.OUT)   
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
        function.addParameter('CVODE_flag', dtype='int32', direction=function.OUT)
        function.addParameter('root_finding_flag', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def a_out_div_a_in_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('m1', dtype='float64', direction=function.IN)
        function.addParameter('m2', dtype='float64', direction=function.IN)
        function.addParameter('m3', dtype='float64', direction=function.IN)
        function.addParameter('e_out', dtype='float64', direction=function.IN)
        function.addParameter('itot', dtype='float64', direction=function.IN)        
        function.result_type = 'float64'
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
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_relative_tolerance():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "Relative tolerance, default 1e-10")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_input_precision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.OUT,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_input_precision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='float64',direction=function.IN,description = "")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_linear_solver():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "linear solver")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_linear_solver():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "linear solver")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_equations_of_motion_specification():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='int32',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_dynamical_stability():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_collision():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_inner_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "...")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_check_for_outer_RLOF():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "...")
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_quadrupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include quadrupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function         

    @legacy_function
    def set_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_octupole_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include octupole terms in the secular equations of motion")
        function.result_type = 'int32'
        return function   

    @legacy_function
    def get_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_1PN_inner_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 1PN 'interaction terms' in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_inner_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 2.5PN inner binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_include_25PN_outer_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Include 2.5PN outer binary terms in the equations of motion")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_inner_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the inner binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_include_outer_tidal_terms():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_ignore_tertiary():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function    
        
    @legacy_function
    def get_ignore_tertiary():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "Whether or not to include the effects of tides in the outer binary system")
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_verbose():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.IN,description = "..")
        function.result_type = 'int32'
        return function    

    @legacy_function
    def get_verbose():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool',direction=function.OUT,description = "..")
        function.result_type = 'int32'
        return function


class SecularTriple(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  SecularTripleInterface(**options), **options)
        self.model_time = 0.0 | units.Myr
        self.evolve_further_after_root_was_found = False ### in some cases, it is desirable to integrate until the given end time, despite the finding of a root
        self.take_into_account_RLOF_after_no_longer_filling_Roche_lobe = True ### this may seem like a strange `feature', but exists for consistency with the stellar/binary evolution code

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
            default_value = 1.0e-10
        )
        object.add_method_parameter(
            "get_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "set_threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "threshold_value_of_e_in_for_setting_tidal_e_in_dot_zero",
            "",
            default_value = 1.0e-12
        )
        object.add_method_parameter(
            "get_input_precision",
            "set_input_precision",
            "input_precision",
            "Input_precision, default 1e-5",
            default_value = 1.0e-5
        )        
        object.add_method_parameter(
            "get_linear_solver",
            "set_linear_solver",
            "linear_solver",
            "linear_solver, default 0",
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
            "get_ignore_tertiary",
            "set_ignore_tertiary",
            "ignore_tertiary",
            "", 
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
            "get_verbose",
            "set_verbose",
            "verbose",
            "..", 
            default_value = True
        )    
    def define_methods(self, object):
        unit_lum = unit_m*unit_l**2/(unit_t**3)
        object.add_method(
            "evolve",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                unit_l,                     ### R1
                unit_l,                     ### R2
                unit_l,                     ### R3
                1.0/unit_t,                 ### spin_angular_frequency1
                1.0/unit_t,                 ### spin_angular_frequency2
                1.0/unit_t,                 ### spin_angular_frequency3                                
                object.NO_UNIT,             ### AMC_star1
                object.NO_UNIT,             ### AMC_star2
                object.NO_UNIT,             ### AMC_star3
                object.NO_UNIT,             ### gyration_radius_star1
                object.NO_UNIT,             ### gyration_radius_star2
                object.NO_UNIT,             ### gyration_radius_star3
                1.0/units.s,                ### k_div_T_tides_star1
                1.0/units.s,                ### k_div_T_tides_star2
                1.0/units.s,                ### k_div_T_tides_star3                                
                unit_l,                     ### a_in
                unit_l,                     ### a_out
                object.NO_UNIT,             ### e_in
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### INCL_in
                object.NO_UNIT,             ### INCL_out                
                object.NO_UNIT,             ### AP_in
                object.NO_UNIT,             ### AP_out
                object.NO_UNIT,             ### LAN_in
                object.NO_UNIT,             ### LAN_out
                unit_t,                     ### t
                unit_t,                     ### dt
            ),
            (
                unit_m,                     ### m1_output
                unit_m,                     ### m2_output
                unit_m,                     ### m3_output
                unit_l,                     ### R1_output
                unit_l,                     ### R2_output
                unit_l,                     ### R3_output
                unit_e/unit_t,              ### tidal_E1_dot_output
                unit_e/unit_t,              ### tidal_E2_dot_output
                unit_e/unit_t,              ### tidal_E3_dot_output
                1.0/unit_t,                 ### spin_angular_frequency1_output
                1.0/unit_t,                 ### spin_angular_frequency2_output
                1.0/unit_t,                 ### spin_angular_frequency3_output                         
                unit_l,                     ### a_in_output
                unit_l,                     ### a_out_output
                object.NO_UNIT,             ### e_in_output
                object.NO_UNIT,             ### e_out_output
                object.NO_UNIT,             ### INCL_in_output
                object.NO_UNIT,             ### INCL_out_output
                object.NO_UNIT,             ### INCL_in_out_output
                object.NO_UNIT,             ### AP_in_output
                object.NO_UNIT,             ### AP_out_output                                
                object.NO_UNIT,             ### LAN_in_output
                object.NO_UNIT,             ### LAN_out_output                                
                unit_t,                     ### t_output
                object.NO_UNIT,             ### CVODE_flag
                object.NO_UNIT,             ### root_finding_flag
                object.ERROR_CODE,
            )
        )

        object.add_method(
            "a_out_div_a_in_dynamical_stability",
            (
                unit_m,                     ### m1
                unit_m,                     ### m2
                unit_m,                     ### m3
                object.NO_UNIT,             ### e_out
                object.NO_UNIT,             ### itot -- NOTE: should be in radians
            ),
            (
                object.NO_UNIT,             ### a_out/a_in for dynamical stability
            )
        )

    def define_particle_sets(self,object):
        object.define_inmemory_set('particles')

    def evolve_model(self,end_time):
        if end_time is None:
            if self.parameters.verbose == True:
                print 'SecularTriple -- please specify end time!'
            return

        parameters = self.parameters
        self.time_step = end_time - self.model_time            

        ### extract data from particles ###
        particles = self.particles
        binaries = particles[particles.is_binary]
        inner_binary = binaries[0]
        outer_binary = binaries[1]
        star1 = inner_binary.child1
        star2 = inner_binary.child2
        if outer_binary.child1.is_binary == True:
            star3 = outer_binary.child2
        elif outer_binary.child2.is_binary == True:
            star3 = outer_binary.child1
            
        m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out = give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary)

        ### the following quantities are not used if tides are not taken into account, but need to pass them to evolve() ###
        spin_angular_frequency1=spin_angular_frequency2=spin_angular_frequency3 = 0.0 | 1.0/units.s
        AMC_star1 = AMC_star2 = AMC_star3 = 0.0
        gyration_radius_star1 = gyration_radius_star2 = gyration_radius_star3 = 0.0
        k_div_T_tides_star1 = k_div_T_tides_star2 = k_div_T_tides_star3 = 0.0 | 1.0/units.s

        if parameters.include_inner_tidal_terms == True:
            try:
                stellar_type1 = star1.stellar_type
                stellar_type2 = star2.stellar_type

                luminosity_star1 = star1.luminosity
                luminosity_star2 = star2.luminosity

                m1_envelope = star1.convective_envelope_mass
                m2_envelope = star2.convective_envelope_mass

                R1_envelope = star1.convective_envelope_radius
                R2_envelope = star2.convective_envelope_radius
                
                AMC_star1 = star1.apsidal_motion_constant
                AMC_star2 = star2.apsidal_motion_constant
                gyration_radius_star1 = star1.gyration_radius
                gyration_radius_star2 = star2.gyration_radius

                spin_angular_frequency1 = star1.spin_angular_frequency
                spin_angular_frequency2 = star2.spin_angular_frequency
                
                k_div_T_tides_star1 = tidal_friction_constant.tidal_friction_constant(stellar_type1,m1,m2,a_in,R1,m1_envelope,R1_envelope,luminosity_star1,spin_angular_frequency1,gyration_radius_star1)
                k_div_T_tides_star2 = tidal_friction_constant.tidal_friction_constant(stellar_type2,m2,m1,a_in,R2,m2_envelope,R2_envelope,luminosity_star2,spin_angular_frequency2,gyration_radius_star2)
        
            except AttributeError:
                print "More attributes required for inner tides"
                exit(-1)
                
        if parameters.include_outer_tidal_terms == True:
            try:
                stellar_type3 = star3.stellar_type

                luminosity_star3 = star3.luminosity

                m3_envelope = star3.convective_envelope_mass

                R3_envelope = star3.convective_envelope_radius
                
                AMC_star3 = star3.apsidal_motion_constant
                gyration_radius_star3 = star3.gyration_radius

                spin_angular_frequency3 = star3.spin_angular_frequency
                
                k_div_T_tides_star3 = tidal_friction_constant.tidal_friction_constant(stellar_type1,m3,m1+m2,a_out,R3,m3_envelope,R3_envelope,luminosity_star3,spin_angular_frequency3,gyration_radius_star3)
        
            except AttributeError:
                print "More attributes required for outer tides"
                exit(-1)
                
        args = (m1,m2,m3,R1,R2,R3, \
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3, \
            AMC_star1,AMC_star2,AMC_star3, \
            gyration_radius_star1,gyration_radius_star2,gyration_radius_star3, \
            k_div_T_tides_star1,k_div_T_tides_star2,k_div_T_tides_star3, \
            a_in,a_out, \
            e_in,e_out, \
            INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out, \
            self.model_time,self.time_step)
            
        ### solve system of ODEs ###
        m1,m2,m3,R1,R2,R3, \
            tidal_E1_dot,tidal_E2_dot,tidal_E3_dot, \
            spin_angular_frequency1,spin_angular_frequency2,spin_angular_frequency3, \
            a_in,a_out,e_in,e_out, \
            INCL_in,INCL_out,INCL_in_out, \
            AP_in,AP_out,LAN_in,LAN_out, \
            end_time_cvode,CVODE_flag,root_finding_flag = self.evolve(*args)
        if self.parameters.verbose == True:
            print 'SecularTriple -- done; t/Myr = ', end_time.value_in(units.Myr),' a_in/AU=',a_in.value_in(units.AU),' a_out/AU=',a_out.value_in(units.AU),'; e_in=',e_in,'e_out=',e_out,' rel_INCL = ',INCL_in,' spin_freq1/day^{-1} = ',spin_angular_frequency1.value_in(1.0/units.day),'LAN_in',LAN_in,'LAN_out',LAN_out,'AP_in',AP_in,'AP_out',AP_out

        print_CVODE_output(self,CVODE_flag)
        self.flag = CVODE_flag
        self.root_finding_flag = root_finding_flag

        ####################
        ### root finding ###
        ####################
        self.dynamical_instability = False
        self.inner_collision = False
        self.outer_collision = False
        self.RLOF_star1 = False
        self.RLOF_star2 = False
        self.RLOF_star3 = False
        
        if CVODE_flag==CV_ROOT_RETURN:
            
            ### dynamical instability ###
            if root_finding_flag==1:
                if self.parameters.verbose == True:
                    print 'SecularTriple -- triple dynamical instability'
                self.dynamical_instability = True

            ### inner collision ###
            if root_finding_flag==2: 
                self.inner_collision = True
                if self.parameters.verbose == True:
                    print 'SecularTriple --inner collision'
            ### outer collision ###
            if root_finding_flag==3: 
                self.outer_collision = True
                if self.parameters.verbose == True:
                    print 'SecularTriple --outer collision'

            ### RLOF star1 ###
            if root_finding_flag==4: 
                self.RLOF_star1 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 1 has filled its Roche Lobe during secular integration'
            if root_finding_flag==-4: 
                self.RLOF_star1 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 1 no longer fills its Roche Lobe during secular integration'

            ### RLOF star2 ###
            if root_finding_flag==5: 
                self.RLOF_star2 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 2 has filled its Roche Lobe during secular integration'
            if root_finding_flag==-5: 
                self.RLOF_star2 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 2 no longer fills its Roche Lobe during secular integration'

            ### RLOF star3 ###
            if root_finding_flag==6: 
                self.RLOF_star3 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 3 has filled its Roche Lobe during secular integration'
            if root_finding_flag==-6: 
                self.RLOF_star3 = True
                if self.parameters.verbose == True:
                    print 'SecularTriple -- star 3 no longer fills its Roche Lobe during secular integration'


        ### update model time ###
        self.model_time += end_time_cvode

        ### update particles ###
        if parameters.include_inner_tidal_terms == True:
            star1.spin_angular_frequency = spin_angular_frequency1
            star2.spin_angular_frequency = spin_angular_frequency2

            star1.tidal_E_dot = tidal_E1_dot
            star2.tidal_E_dot = tidal_E2_dot

        if parameters.include_outer_tidal_terms == True:
            star3.spin_angular_frequency = spin_angular_frequency3
            
            star3.tidal_E_dot = tidal_E3_dot

        inner_binary.semimajor_axis = a_in
        inner_binary.eccentricity = e_in
        outer_binary.semimajor_axis = a_out
        outer_binary.eccentricity = e_out
        inner_binary.inclination = INCL_in
        outer_binary.inclination = INCL_out
        inner_binary.argument_of_pericenter = AP_in
        outer_binary.argument_of_pericenter = AP_out
        inner_binary.longitude_of_ascending_node = LAN_in
        outer_binary.longitude_of_ascending_node = LAN_out


def print_CVODE_output(self,CVODE_flag):
    if self.parameters.verbose == False: return
    
    if CVODE_flag==CV_SUCCESS:
        print "SecularTriple -- ODE integration proceeded successfully."
    elif CVODE_flag==CV_ROOT_RETURN:
        print "SecularTriple -- root was found during ODE integration."
    elif CVODE_flag==CV_WARNING:
        print "SecularTriple -- (recoverable) warnings occurred during ODE integration."
    if CVODE_flag<0:
        if CVODE_flag==CV_TOO_MUCH_WORK:
            message = "too many steps were taken during ODE integration."
        elif CVODE_flag==CV_TOO_MUCH_ACC:
            message = "required accuracy could not be obtained."
        elif CVODE_flag==CV_ERR_FAILURE:
            message = "error test failures occurred too many times during one internal time-step or minimum step size was reached."
        elif CVODE_flag==CV_CONV_FAILURE:
            message = "convergence test failures occurred too many times during one internal time-step or minimum step size was reached."
        elif CVODE_flag==CV_LINIT_FAIL:    
            message = "the linear solver's initialization function failed."
        elif CVODE_flag==CV_LSETUP_FAIL:    
            message = "the linear solver's setup function failed in an unrecoverable manner."
        elif CVODE_flag==CV_LSOLVE_FAIL:    
            message = "the linear solver's solve function failed in an unrecoverable manner."
        elif CVODE_flag==CV_RHSFUNC_FAIL:
            message = "the right-hand side function failed in an unrecoverable manner."
        elif CVODE_flag==CV_FIRST_RHSFUNC_ERR:    
            message = "the right-hand side function failed at the first call."
        elif CVODE_flag==CV_REPTD_RHSFUNC_ERR:    
            message = "the right-hand side function had repetead recoverable errors."
        elif CVODE_flag==CV_UNREC_RHSFUNC_ERR:    
            message = "the right-hand side function had a recoverable error, but no recovery is possible."
        elif CVODE_flag==CV_RTFUNC_FAIL:    
            message = "the rootfinding function failed in an unrecoverable manner."
        elif CVODE_flag==CV_MEM_FAIL:    
            message = "a memory allocation failed."
        elif CVODE_flag==CV_MEM_NULL:    
            message = "the cvode mem argument was NULL."
        elif CVODE_flag==CV_ILL_INPUT:    
            message = "one of the function inputs is illegal."
        elif CVODE_flag==CV_NO_MALLOC:    
            message = "the cvode memory block was not allocated by a call to CVodeMalloc."
        elif CVODE_flag==CV_BAD_K:    
            message = "the derivative order k is larger than the order used."
        elif CVODE_flag==CV_BAD_T:    
            message = "the time t is outside the last step taken."
        elif CVODE_flag==CV_BAD_DKY:    
            message = "the output derivative vector is NULL."
        elif CVODE_flag==CV_TOO_CLOSE:    
            message = "the output and initial times are too close to each other."
        else:
            message = "unknown error."
        print "SecularTriple -- unrecoverable error occurred during secular integration (CVODE_flag ",str(CVODE_flag),"): ",message

def give_stellar_masses_radii_and_binary_parameters(self,star1,star2,star3,inner_binary,outer_binary):
    
    ### stellar parameters ###
    m1 = star1.mass
    m2 = star2.mass
    m3 = star3.mass
    
    R1 = star1.radius
    R2 = star2.radius
    R3 = star3.radius

    ### binary parameters
    a_in = inner_binary.semimajor_axis
    a_out = outer_binary.semimajor_axis
    e_in = inner_binary.eccentricity
    e_out = outer_binary.eccentricity

    INCL_in = inner_binary.inclination
    INCL_out = outer_binary.inclination
    
    AP_in = inner_binary.argument_of_pericenter
    AP_out = outer_binary.argument_of_pericenter
    LAN_in = inner_binary.longitude_of_ascending_node
    LAN_out = outer_binary.longitude_of_ascending_node
    
    return m1,m2,m3,R1,R2,R3,a_in,a_out,e_in,e_out,INCL_in,INCL_out,AP_in,AP_out,LAN_in,LAN_out

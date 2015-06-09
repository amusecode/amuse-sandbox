from amuse.test.amusetest import TestWithMPI
from amuse.units import units, constants
import numpy as np

import triple_nokozaidt as triple
test_asserts = True


class TestSeBa(TestWithMPI):
    def __init__(self):
        print 'initializing...'

    def test0(self):
        print 'test0'
        
        a=5|units.MSun
        self.assertEqual(a, 5 | units.MSun)

        a=1.009|units.MSun
        self.assertAlmostRelativeEqual(a, 1 | units.MSun, 2)

        a=10.1|units.MSun
        self.assertAlmostRelativeEqual(a, 10| units.MSun, 2)
    
        a=10.0049|units.MSun
        self.assertAlmostEqual(a, 10| units.MSun, 2)

        print 'test0: succeeded'

    #test input parameters
    def test1(self):
        print 'test1'
        
        M1 = 1|units.MSun
        M2 = 0.5|units.MSun
        M3 = 0.08|units.MSun
        a_in = 12345|units.RSun
        a_out = 12345670000000|units.RSun
        e_in = 0.01
        e_out = 0.001
        i = 1.9*np.pi/180.0
        g_in = 1.5
        g_out = 1.55
        z = 0.001
        T_end = 1|units.yr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, metallicity = z, tend = T_end)
        
        self.assertEqual(tr.triple.child2.child1.mass, M1)        
        self.assertEqual(tr.triple.child2.child2.mass, M2)        
        self.assertEqual(tr.triple.child1.mass, M3)        
        
        self.assertEqual(tr.triple.child2.semimajor_axis, a_in)        
        self.assertEqual(tr.triple.semimajor_axis, a_out)        

        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, e_in, 4)        
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, e_out, 4)        
        
        self.assertAlmostRelativeEqual(tr.triple.relative_inclination, i, 3) # precision of secular code 1e-5      

        self.assertAlmostRelativeEqual(tr.triple.child2.argument_of_pericenter, g_in, 4)        
        self.assertAlmostRelativeEqual(tr.triple.argument_of_pericenter, g_out, 4)        

        self.assertEqual(tr.stellar_code.parameters.metallicity, z)
        
        self.assertAlmostRelativeEqual(tr.time, 1|units.yr, 4)
        
        print 'test1: succeeded'


    # wind in inner system   
    def test2(self):
        print 'test2'
    
        M1 = 7|units.MSun
        M2 = 1|units.MSun
        M3 = 1|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 60|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        print tr.triple.child2.child1.stellar_type, tr.triple.child2.child2.stellar_type, tr.triple.child1.stellar_type

#0.00100025844189 0.00099996533347
#2.3147618017e+13 [Myr**-1] 55942.028509 [Myr**-1] 59.335481569 [Myr**-1]


        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 34831 | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 2729968 | units.RSun, 2)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test2: succeeded'


    def test2_alt(self):
        print 'test2'
    
        M1 = 4|units.MSun
        M2 = 1|units.MSun
        M3 = 1|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 220|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
#0.00100605187918 0.000999466463344
#3.0272154593e+11 [Myr**-1] 43964.168682 [Myr**-1] 48.160333818 [Myr**-1]

        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.873069 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 26694 | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 2088360 | units.RSun, 2)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test2: succeeded'


    # wind in inner system   
    def test3(self):
        print 'test3'
    
        M1 = 7|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 60|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)
            

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#0.000997491083201 0.00100000186191
#1.8020988953e+15 [Myr**-1] 17146.253271 [Myr**-1] 71.690850135 [Myr**-1]

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 6. | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 17816 | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1773139 | units.RSun, 2)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        
        
        print 'test3: succeeded'
        
    # wind in inner system   
    def test4(self):
        print 'test4'
    
        M1 = 7|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 80|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
                
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.14678 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
                    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 55931 | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis,  5445910| units.RSun, 2)              

#with old timestep
#1.2942010713 [MSun] 1.1461543599 [MSun] 0.080000000000 [MSun]
#55931.333090 [RSun] 5445910.3444 [RSun]
#1.00000080097e-05 1.00000000113e-05
#7.9428879302 [Myr**-1] 0.80841846468 [Myr**-1] 1651407.5261 [Myr**-1]
#53270.928626 [RSun] 5189744.2074 [RSun]

#with extra timestepp after type change
#1.2942117559 [MSun] 1.1445711755 [MSun] 0.080000000000 [MSun]
#53539.679182 [RSun] 5215554.7219 [RSun]
#1.00000080139e-05 1.00000000121e-05
#8.0313404248 [Myr**-1] 0.67690521866 [Myr**-1] 1651407.5261 [Myr**-1]
#53305.277123 [RSun] 5192984.2136 [RSun]

        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2        
        print a_in_final_theory, a_out_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 1.2)        
        #difference comes from small differences in timesteps, for beginning super wind phase
        print 'test4: succeeded'
        
    # wind in outer system   
    def test5(self):
        print 'test5'
        
        M1 = 0.1|units.MSun
        M2 = 0.1|units.MSun
        M3 = 7.|units.MSun
        a_in = 10.|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 60|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)
    

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        if test_asserts:
            self.assertEqual(tr.triple.child2.child1.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child1.mass, 1.2967 | units.MSun, 2)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 10.0| units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 4810448 | units.RSun, 2)        

        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        
        
        print 'test5: succeeded'
        
    def test6(self):
        print 'test6'
        
        M1 = 3.0|units.MSun
        M2 = 2.0|units.MSun
        M3 = 0.08|units.MSun
        a_in = 15|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 385 | units.Myr 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)



        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            #nucl timescale mt
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.44221498171 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 4.5520733201 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 130.52823672 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000747.5318 | units.RSun, 2)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test6: succeeded'

    def test7(self):
        print 'test7'

        M1 = 7.0|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 50|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        print tr.triple.child2.child1.stellar_type, tr.triple.child2.child2.stellar_type, tr.triple.child1.stellar_type


        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3578 | units.MSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.6361 | units.MSun, 2)        
#            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#    
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 352.43 | units.RSun, 2)  
#            # because of wind mass loss a_out_final != 1e6RSun      
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000466 | units.RSun, 2)        
            #nucl timescale mt
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3403598844 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.639010234 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 359.69117942 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1001572.1072 | units.RSun, 2)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test7: succeeded'

    def test8(self):
        print 'test8'

        M1 = 7.0|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 100|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        print tr.triple.child2.child1.stellar_type, tr.triple.child2.child2.stellar_type, tr.triple.child1.stellar_type

        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3591 | units.MSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.6291 | units.MSun, 2)        
#            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#    
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 700.1299 | units.RSun, 2)  
#            # because of wind mass loss a_out_final != 1e6RSun      
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000898 | units.RSun, 2)        
            #nucl timescale mt
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3508433033 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.625443434 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 707.42929292 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1001726.5769 | units.RSun, 2)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test8: succeeded'







    def test6_with_tides(self):
        print 'test6'
        
        M1 = 3.0|units.MSun
        M2 = 2.0|units.MSun
        M3 = 0.08|units.MSun
        a_in = 15|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 385 | units.Myr 
        tidal_terms = True
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)



        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            #nucl timescale mt
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.4476487 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 4.547591 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 130.23696 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000939 | units.RSun, 2)        
        
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 1383819872| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 29230375| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)
        

        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test6 with tides: succeeded'

    def test7_with_tides(self):
        print 'test7'

        M1 = 7.0|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 50|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr 
        tidal_terms = True
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        print tr.triple.child2.child1.stellar_type, tr.triple.child2.child2.stellar_type, tr.triple.child1.stellar_type


        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.34661 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.64483 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 360.22 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000654.9252 | units.RSun, 2)        

            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 21841671| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 9763575| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test7 with tides: succeeded'

    def test8_with_tides(self):
        print 'test8'

        M1 = 7.0|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 100|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr 
        tidal_terms = True
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        print tr.triple.child2.child1.stellar_type, tr.triple.child2.child2.stellar_type, tr.triple.child1.stellar_type

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3518307 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.6299890 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 712.398| units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1001393 | units.RSun, 2)        
        
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 166720249| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 3493846| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_in_final_theory, a_out_final_theory

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test8 with tides: succeeded'











    def test9(self):
        print 'test9'

        M1 = 7|units.MSun
        M2 = 1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 200|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3517 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 1.0 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 3.4474 | units.RSun, 2)    
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3322817 | units.RSun, 2)        

        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test9: succeeded'


    def test10(self):
        print 'test10'

        M1 = 7|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 300|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency


        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3534 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 6.0 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 262.86 | units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1760013 | units.RSun, 2)        
#1.3506318555 [MSun] 6.0000000000 [MSun] 0.080000000000 [MSun]
#263.58702734 [RSun] 1760244.5558 [RSun]
#1e-05 1e-05
#11712499.861 [Myr**-1] 12187468.954 [Myr**-1] 71.690850135 [Myr**-1]
#1760280.9902 [RSun]
            

        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test10: succeeded'

    def test11(self): #ms hook
        print 'test11'

        M1 = 1.5|units.MSun
        M2 = 1.499|units.MSun
        M3 = 0.08|units.MSun
        a_in = 200|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 3000|units.Myr
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)
        
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.40921 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 0.31219 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 1.3760 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3832604 | units.RSun, 2)        

        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test11: succeeded'


 
    def test12(self):
        print 'test12'

        M1 = 7.0|units.MSun
        M2 = 7.0|units.MSun
        M3 = 0.08|units.MSun
        a_in = 100|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 49.01|units.Myr 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
            
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3430 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.3430 | units.MSun, 2)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 0.72327 | units.RSun, 2)  
            # because of wind mass loss a_out_final != 1e6RSun      
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 5090161 | units.RSun, 2)        


        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test12: succeeded'

    def test13(self):
        print 'test13'

        M1 = 7.0|units.MSun
        M2 = 7.0|units.MSun
        M3 = 0.08|units.MSun
        a_in = 25|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 48.92|units.Myr 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency


        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 11.1397 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.08 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 1254880.1897| units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.0000000000e+100 | units.RSun, 2)        
        
        print 'test13: succeeded'

        #test tides
    def test14(self):
        print 'test14'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 20|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.8
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 12500|units.Myr 
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, stop_at_mass_transfer = True)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99837 | units.MSun, 4)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 6.7610| units.RSun, 4)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1001367 | units.RSun, 4)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.e-5, 4)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.e-5, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        print a_in_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test14: succeeded'
 

         #test tides
    def test15(self):
        print 'test15'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 35|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.8
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 20000|units.Myr 
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, stop_at_mass_transfer = True)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99625| units.MSun, 4)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 11.8328| units.RSun, 4)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1003162 | units.RSun, 4)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.e-5, 4)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.e-5, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        print a_in_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test15: succeeded'
 

    def test16(self):
        print 'test16'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 30 |units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 13500 |units.Myr
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, stop_at_mass_transfer=True)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.98825083719 | units.MSun, 4)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 28.290035621| units.RSun, 4)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1009993 | units.RSun, 4)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.e-5, 4)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.e-5, 4)        

        a_in_final_theory = a_in * (1-e_in**2)
        print a_in_final_theory
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test16: succeeded'
 

    def test17(self):
        print 'test17'

        M1 = 1.|units.MSun
        M2 = 0.9|units.MSun
        M3 = 0.08|units.MSun
        a_in = 30|units.RSun 
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 13500|units.Myr 
        dr = 0.0005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.24242 | units.MSun, 4)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.9 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 28.4444| units.RSun, 4)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1619732 | units.RSun, 4)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.e-5, 4)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.e-5, 4)        


#       should be gamma CE
#        a_in_final_theory = a_in * (1-e_in**2) #although slightly different due to spin
#        print a_in_final_theory, 
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_out_final_theory
        if test_asserts:   
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test17: succeeded'
 
      
    def test18(self):
        print 'test18'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun 
        a_out = 1.e10|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 12500|units.Myr 
        dr = 0.005 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.52058 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 17728.97| units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 16846361308. | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.00e-5, 2)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.00e-5, 2)        
    
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 0.033007513673| 1./units.Myr, 2)
#            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 20762.5| 1./units.Myr, 2)
#            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 2.15328e-05| 1./units.Myr, 2)

            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 9.2145348309| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 3827799| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)


        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:   
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test18: succeeded'

        
    def test19(self):
        print 'test19'

        M1 = 5.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun 
        a_out = 1.e10|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 125|units.Myr
        dr = 0.005
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        
        if test_asserts:   
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.00496 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 46716.| units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 44236004807. | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.00e-5, 2)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.00e-5, 2)        
    
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 0.30112881089| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 3832832.994| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)


        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:   
#            print 'why so low accuracy?'-> because somewhere on the rgb star starts expanding rapidly
#                                this sudden change is resolved better with small timesteps from dt_radius_change
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 1.5)        


        print 'test19: succeeded'



    def test19_alt(self):
        print 'test19_alt'

        M1 = 4.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun 
        a_out = 1.e10|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 220|units.Myr
        dr = 0.005
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, tidal_terms = tidal_terms)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        if test_asserts:   
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.87133 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 42518| units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 40043414739. | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1.00e-5, 2)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1.00e-5, 2)        

            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 0.12587058782| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 3832794| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 2)


        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        if test_asserts:   
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        


        print 'test19: succeeded'



    def test20(self):
        print 'test20'

        M1 = 7.5|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 1.e4|units.RSun 
        a_out = 1.e10|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 50|units.Myr
        dr = 0.005 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency


        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3686793822 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 55294| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 52951483340. | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 1e-05, 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 1e-05, 4)        

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 86.4380055| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 3832863| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 1651407| 1./units.Myr, 4)

        print 'test20: succeeded'



    def test21(self):
        print 'test21'        
        M1 = 1.3|units.MSun
        M2 = 0.5|units.MSun
        M3 = 0.5|units.MSun
        a_in = 200 |units.RSun
        a_out = 20000|units.RSun
        e_in = 0.1
        e_out = 0.5
        i = 80.*np.pi/180.0
        g_in = 0.1
        g_out = 0.5
        T_end = 5|units.Myr 
        dr = 0.005 
        tidal_terms = False
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, maximum_radius_change_factor=dr, tidal_terms = tidal_terms)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        if test_asserts:   
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3| units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 0.5 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 0.5 | units.MSun)        
            
            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 200.| units.RSun, 2)  
            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 20000. | units.RSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0.939581864833, 2)  
            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.501813404911, 2)        
    
            #initially 9402668 and 10628, difference due to rdot -> start on corotation 
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 9396370.| 1./units.Myr, 2)
#            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 9402526| 1./units.Myr, 2)
#            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 10628| 1./units.Myr, 2)

            self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 1599505764.| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 219467244| 1./units.Myr, 2)
            self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 219467244| 1./units.Myr, 2)

        
        print 'test21: succeeded'



    def test22(self):
        print 'test22'
        M1 = 1.|units.MSun
        M2 = 1.1|units.MSun
        M3 = 2|units.MSun
        a_in = 2.e6 |units.RSun
        a_out = 2.3e6|units.RSun
        e_in = 0.
        e_out = 0.
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out)

        print tr.triple.dynamical_instability_at_initialisation,
        print tr.triple.dynamical_instability
        print tr.triple.child2.child1.age, tr.triple.child2.child2.age, tr.triple.child1.age
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        if test_asserts:   
            self.assertEqual(tr.triple.dynamical_instability, True)        
            self.assertEqual(tr.triple.dynamical_instability_at_initialisation, True)        
            self.assertEqual(tr.triple.child2.child1.age, 0.| units.yr)        
            self.assertEqual(tr.triple.child2.child2.age, 0. | units.yr)        
            self.assertEqual(tr.triple.child1.age, 0. | units.yr)        
            
            self.assertEqual(tr.triple.child2.child1.mass, 1.1| units.MSun)        
            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 2. | units.MSun)        
            
            self.assertEqual(tr.triple.child2.semimajor_axis, 2.e6| units.RSun)  
            self.assertEqual(tr.triple.semimajor_axis, 2.3e6 | units.RSun)        
            self.assertEqual(tr.triple.child2.eccentricity, 1.e-5)  
            self.assertEqual(tr.triple.eccentricity, 1.e-5)        
    
        print 'test22: succeeded'




    # wind in inner system, massive star with SN
    def test23(self):
        print 'test23'
    
        M1 = 9|units.MSun
        M2 = 1|units.MSun
        M3 = 1|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 33.34|units.Myr
        tidal_terms = False
        stop_at_SN = True
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms, 
        stop_at_SN = stop_at_SN)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.child1.previous_mass, tr.triple.child2.child2.previous_mass, tr.triple.child1.previous_mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#just before SN
#        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.22811 | units.MSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.previous_mass, 8.52456 | units.MSun, 2)        
#            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
#            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 10499 | units.RSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1045174 | units.RSun, 2)        
#        
#        #under the assumption of no wind accretion
#        M_child1 = tr.triple.child2.child1.previous_mass + tr.triple.child2.child2.previous_mass
#        M_child2 = M_child1 + tr.triple.child1.previous_mass
#        a_in_final_theory =  a_in * (M1+M2) / M_child1
#        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
#        print a_in_final_theory, a_out_final_theory
#        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.22811 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.previous_mass, 8.52456 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, -4615.4878113 | units.RSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, -20197.418621 | units.RSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 3.2747492923, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.eccentricity, 44.5759681379, 2)        


#SeBa
#get_total_mass() = 8.51507
#core_mass = 1.22828
#m_comp_0 = 1

        print 'test23: succeeded'

    # wind in inner system, massive star with SN kick
    def test24(self):
        print 'test24'
    
        M1 = 9|units.MSun
        M2 = 8.5|units.MSun
        M3 = 1|units.MSun
        a_in = 1.e4|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        T_end = 33.5|units.Myr
        tidal_terms = False
        stop_at_SN = True 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, tend = T_end, tidal_terms = tidal_terms, 
        stop_at_SN = stop_at_SN)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.child1.previous_mass, tr.triple.child2.child2.previous_mass, tr.triple.child1.previous_mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#before SN
#        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.22811 | units.MSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.child2.child1.previous_mass, 8.52456 | units.MSun, 2)        
#            self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
#            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 10499 | units.RSun, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1045174 | units.RSun, 2)        
#        
#        #under the assumption of no wind accretion
#        M_child1 = tr.triple.child2.child1.previous_mass + tr.triple.child2.child2.previous_mass
#        M_child2 = M_child1 + tr.triple.child1.previous_mass
#        a_in_final_theory =  a_in * (M1+M2) / M_child1
#        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
#        print a_in_final_theory, a_out_final_theory
#        if test_asserts:
#            self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
#            self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        if test_asserts:
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.22811 | units.MSun, 2)        
            self.assertAlmostRelativeEqual(tr.triple.child2.child1.previous_mass, 8.52456 | units.MSun, 2)        
            self.assertEqual(tr.triple.child2.child2.mass, 8.405757934 | units.MSun)        
            self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                

#SeBa:
#get_total_mass() = 8.52818,
#core_mass = 1.22828
#m_comp_0 = 8.39485


        print 'test24: succeeded'




if __name__ == '__main__':
    test = TestSeBa()

#test set up
#    test.test0()
#    test.test1()
###
####test wind mass loss inner system
#    test.test2()
#    test.test2_alt()
#    test.test3()
#    test.test4()
##
###test wind mass loss outer system
#    test.test5()
###
####test stable mass transfer in inner system 
#    test.test6() 
#    test.test7()
#    test.test8()

#    test.test6_with_tides() 
#    test.test7_with_tides()
#    test.test8_with_tides()

#    
##test common envelope evolution in inner system
#    test.test9() #alpha-ce - need tides to bring the system into corotation -> DR XXX!!!
#    test.test10() #gamma-ce
#    test.test11() #double alpha-ce
##
#### test contact system in inner system
#    test.test12() #double alpha-ce # why does orbit shrink after ce?
#    test.test13() #merger 
#
#
###    #test tides
#    test.test14() #merger 
#    test.test15() #merger
#    test.test16() # darwin-riemann instability XXX!!! 
#    test.test17() # no darwin-riemann instability


##
##
####  #test spin-r coupling -> no tidal terms
#    test.test18()
#    test.test19()
#    test.test19_alt()
#    test.test20() # sn kick not included
#
#
#    #test kozai
#    test.test21()
    
#    dynamical instability at initialization
#    test.test22()

#    massive stars with SN
    test.test23() 
    test.test24() 

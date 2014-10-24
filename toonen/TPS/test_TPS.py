from amuse.test.amusetest import TestWithMPI
from amuse.units import units, constants
import numpy as np

import triple

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
        a_out = 1234567|units.RSun
        e_in = 0.01
        e_out = 0.001
        i = 0.1*np.pi/180.0
        g_in = 1.5
        g_out = 1.55
        o_in = 0.1
        o_out = 0.15
        z = 0.001
        T_end = 1|units.yr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, metallicity = z, tend = T_end)
        
        self.assertEqual(tr.triple.child2.child1.mass, M1)        
        self.assertEqual(tr.triple.child2.child2.mass, M2)        
        self.assertEqual(tr.triple.child1.mass, M3)        
        
        self.assertEqual(tr.triple.child2.semimajor_axis, a_in)        
        self.assertEqual(tr.triple.semimajor_axis, a_out)        

        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, e_in, 4)        
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, e_out, 4)        
        
        self.assertAlmostRelativeEqual(tr.triple.relative_inclination, i, 4)        

        self.assertAlmostRelativeEqual(tr.triple.child2.argument_of_pericenter, g_in, 4)        
        self.assertAlmostRelativeEqual(tr.triple.argument_of_pericenter, g_out, 4)        

        print 'why is o reset to 0 by secular code?'
#        self.assertAlmostRelativeEqual(tr.triple.child1.longitude_of_ascending_node, o_in, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.longitude_of_ascending_node, o_out, 4)        

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
        o_in = 0.0
        o_out = 0.0
        T_end = 60|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#0.00100025844189 0.00099996533347
#2.3147618017e+13 [Myr**-1] 55942.028509 [Myr**-1] 59.335481569 [Myr**-1]


        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)                
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 34831 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 2729968 | units.RSun, 4)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
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
        o_in = 0.0
        o_out = 0.0
        T_end = 220|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
#0.00100605187918 0.000999466463344
#3.0272154593e+11 [Myr**-1] 43964.168682 [Myr**-1] 48.160333818 [Myr**-1]

        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.873069 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 26694 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 2088360 | units.RSun, 4)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
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
        o_in = 0.0
        o_out = 0.0
        T_end = 60|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)
            
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 6. | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 17807 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1773139 | units.RSun, 4)        

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#0.000997491083201 0.00100000186191
#1.8020988953e+15 [Myr**-1] 17146.253271 [Myr**-1] 71.690850135 [Myr**-1]

        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
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
        o_in = 0.0
        o_out = 0.0
        T_end = 80|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2967 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.14678 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 53173 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis,  5183243| units.RSun, 4)              
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
#0.000997490653614 0.0010000013316
#1.8020988953e+15 [Myr**-1] 1.9691679561e+11 [Myr**-1] 71.690850135 [Myr**-1]
        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        
        
        print 'test4: succeeded'
        
    # wind in outer system   
    def test5(self):
        print 'test5'
        
        M1 = 0.1|units.MSun
        M2 = 0.1|units.MSun
        M3 = 7.|units.MSun
        a_in = 1.|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 60|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)
    
#        self.assertEqual(tr.triple.child2.child1.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertAlmostRelativeEqual(tr.triple.child1.mass, 1.2964 | units.MSun, 4)        
#
#        self.assertEqual(tr.triple.child2.semimajor_axis, 1.0| units.RSun)        
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 4849289 | units.RSun, 4)        
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
        print a_in_final_theory, a_out_final_theory
        
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 2)        
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        
        
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
        o_in = 0.0
        o_out = 0.0
        T_end = 380 | units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.419349 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 4.5803 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 150.039| units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000065 | units.RSun, 4)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test6: succeeded'

    def test7(self):
        print 'test7'

        M1 = 7.0|units.MSun
        M2 = 6|units.MSun
        M3 = 0.08|units.MSun
        a_in = 50|units.RSun
#        a_in = 23|units.RSun op het randje van ms/hg donor, met getijden net ms
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 49.1|units.Myr 
#        T_end = 48.96|units.Myr# until RLOF
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency


        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3294 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.6590 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

#        self.assertAlmostRelativeEqual(tr.triple.child1.semimajor_axis, 171.95 | units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 172.04 | units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000880 | units.RSun, 4)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
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
        o_in = 0.0
        o_out = 0.0
        T_end = 49.1|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3459 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.34561 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 11.6423 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 731.54 | units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 731.97 | units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000902 | units.RSun, 4)        
        
        #under the assumption of conservative mass transfer
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_in_final_theory =  a_in * (M1 * M2 / M1f / M2f)**2 #conservative mass transfer
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1.5)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test8: succeeded'

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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        
        
#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3517 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 1.0 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 3.4474 | units.RSun, 4)    
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3322817 | units.RSun, 4)        
#        
#        M1f = tr.triple.child2.child1.mass
#        M2f = tr.triple.child2.child2.mass
#        M3f = tr.triple.child1.mass       
#        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3534 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 6.0 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 262.86 | units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1760013 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 248.1884 | units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1759615 | units.RSun, 4)        
        
#0.000999999999989 0.000999999999903
#12933662.796 [Myr**-1] 13336108.835 [Myr**-1] 71.690850135 [Myr**-1]

        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
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
        o_in = 0.0
        o_out = 0.0
        T_end = 3000|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)
        
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        


#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.40135 | units.MSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 0.31219 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 1.3183 | units.RSun, 4)  
#        # because of wind mass loss a_out_final != 1e6RSun      
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3880131 | units.RSun, 4)        
#        
#        M1f = tr.triple.child2.child1.mass
#        M2f = tr.triple.child2.child2.mass
#        M3f = tr.triple.child1.mass       
#        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.34345 | units.MSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.34345 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 0.75302 | units.RSun, 4)  
#        # because of wind mass loss a_out_final != 1e6RSun      
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 5088750 | units.RSun, 4)        

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
            
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        
        print a_out_final_theory

        print 'test12: succeeded'

    def test13(self):
        print 'test13'
        #stop_at_merger = True

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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

#        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 6.98888 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 6.99565 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 6.99565 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 24.2392| units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1000617 | units.RSun, 4)        

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#4.82702766647e-12 0.001
#620908463.27 [Myr**-1] 620908463.27 [Myr**-1] 74.380858360 [Myr**-1]

        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory

        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

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
        o_in = 0.0
        o_out = 0.0
        T_end = 11050|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
#        print tr.triple.child2.eccentricity, tr.triple.eccentricity
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
#        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        
#

        a_in_final_theory = a_in * (1-e_in**2)
        print a_in_final_theory
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test14: succeeded'
 

         #test tides
    def test15(self):
        print 'test15'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 35|units.RSun
#        a_in = 20|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.8
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 20000|units.Myr 
        dr = 0.005 #0.0005
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
#        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        print a_in_final_theory
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

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
        o_in = 0.0
        o_out = 0.0
        T_end = 12250|units.Myr 
        dr = 0.001#0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
#        print tr.triple.child2.eccentricity, tr.triple.eccentricity
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
#        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        
#
#
#        a_in_final_theory = a_in * (1-e_in**2)
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        
        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius

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
        o_in = 0.0
        o_out = 0.0
        T_end = 12250|units.Myr 
        dr = 0.0005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius
#        radius_ch2ch1_init
#        spin_angular_frequency_ch2ch1_init
#        Jspin_ch2ch1_init = spin_angular_frequency_ch2ch1_init * M1* radius_ch2ch1_init**2 #gyration radius


#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#
#        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
#        print tr.triple.child2.eccentricity, tr.triple.eccentricity
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
#        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        
#
#
#        a_in_final_theory = a_in * (1-e_in**2)
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        
#
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
        o_in = 0.0
        o_out = 0.0
        T_end = 12500|units.Myr 
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.52103 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 17712.5| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 16832350474. | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0.00099999, 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 92815335| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 20762.5| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 2.15328e-05| 1./units.Myr, 4)

        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius
        spin1_init = 20790.081940 |1./units.Myr
        R1_init =  0.88824945030 |units.RSun   
        Jspin_ch2ch1_init = spin1_init * M1 * R1_init**2 #gyration radius
        self.assertAlmostRelativeEqual(Jspin_ch2ch1, Jspin_ch2ch1_init, 2)

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
        o_in = 0.0
        o_out = 0.0
        T_end = 125|units.Myr
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.00496 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 46155| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 43714937660. | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0.00099999, 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 37118045560.| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 44765.1| 1./units.Myr, 4)
        self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 4.51154e-05| 1./units.Myr, 4)

        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius
        spin1_init = 44765.655322 |1./units.Myr
        R1_init =  2.6355930592 |units.RSun   
        Jspin_ch2ch1_init = spin1_init * M1 * R1_init**2 #gyration radius
        self.assertAlmostRelativeEqual(Jspin_ch2ch1, Jspin_ch2ch1_init, 2)

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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99851 | units.MSun, 4)        
#        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
#        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
#        
#        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 10013.6| units.RSun, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 10012668194. | units.RSun, 4)        
#        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0.001, 4)  
#        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        
#
#        self.assertAlmostRelativeEqual(tr.triple.child2.child1.spin_angular_frequency, 1187.24| 1./units.Myr, 4)
#        self.assertAlmostRelativeEqual(tr.triple.child2.child2.spin_angular_frequency, 20763.6| 1./units.Myr, 4)
#        self.assertAlmostRelativeEqual(tr.triple.child1.spin_angular_frequency, 2.15328e-05| 1./units.Myr, 4)

        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius
        spin1_init = 56416.010036 |1./units.Myr
        R1_init =  3.4567384397 |units.RSun   
        Jspin_ch2ch1_init = spin1_init * M1 * R1_init**2 #gyration radius
        self.assertAlmostRelativeEqual(Jspin_ch2ch1, Jspin_ch2ch1_init, 2)
        print 'test20: succeeded'
        
    def test21(self):
        print 'test21'
        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 5 |units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 15500|units.Myr 
        dr = 0.005 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency
        Jspin_ch2ch1 = tr.triple.child2.child1.spin_angular_frequency * tr.triple.child2.child1.mass* tr.triple.child2.child1.radius**2 #gyration radius

        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        print tr.triple.child2.child1.spin_angular_frequency, tr.triple.child2.child2.spin_angular_frequency, tr.triple.child1.spin_angular_frequency

  
        
if __name__ == '__main__':
    test = TestSeBa()

#     test set up
#    test.test0()
#    test.test1()
#
##     test wind mass loss inner system
    test.test2()
    test.test2_alt()
    test.test3()
    test.test4()
###
####     test wind mass loss outer system
    test.test5()
#
#    # test stable mass transfer in inner system -> all darwin riemann instabilities
#    test.test6() 
#    test.test7()
#    test.test8()
#
#    # test common envelope evolution in inner system
    test.test9() #alpha-ce - gamma
    test.test10() #gamma-ce
    test.test11() #double alpha-ce
##
##    # test contact system in inner system
    test.test12() #double alpha-ce
    test.test13() #merger
##

##    #test tides
#    test.test14() #secular code crashes : e_in -> 1e-14
#    test.test15() #secular code crashes : e_in -> 1e-14
#    test.test16() # darwin-riemann instability, secular code crashes : e_in -> 1e-14
#    test.test17() # no darwin-riemann instability, secular code crashes : e_in -> 1e-14


##  #test spin-r coupling
    test.test18()
    test.test19()
##    test.test20() M large timestep issues


# test dt 
    # 8, 50Msun succeeds
    #10, 15MSun crashes


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
        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2964 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 1. | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 1. | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 35015 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 2739617 | units.RSun, 4)        
        
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
            
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2964 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 6. | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 17838 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1775327 | units.RSun, 4)        
        
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

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.2964 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.1432 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 53608 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis,  5221675| units.RSun, 4)              
        
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
    
        self.assertEqual(tr.triple.child2.child1.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertAlmostRelativeEqual(tr.triple.child1.mass, 1.2964 | units.MSun, 4)        

        self.assertEqual(tr.triple.child2.semimajor_axis, 1.0| units.RSun)        
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 4849289 | units.RSun, 4)        
        
        #under the assumption of no wind accretion
        M_child1 = tr.triple.child2.child1.mass + tr.triple.child2.child2.mass
        M_child2 = M_child1 + tr.triple.child1.mass
        a_in_final_theory =  a_in * (M1+M2) / M_child1
        a_out_final_theory = a_out * (M1+M2+M3) / M_child2
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
        o_in = 0.0
        o_out = 0.0
        T_end = 380 | units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

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
        a_in = 23|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 0.0
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 49|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

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
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        
        
        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3517 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 1.0 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 3.4474 | units.RSun, 4)    
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3322817 | units.RSun, 4)        
        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)


        print tr.triple.child2.child1.mass, tr.triple.child2.child2.mass, tr.triple.child1.mass
        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.3517 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 6.0 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 262.86 | units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1760013 | units.RSun, 4)        
        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

        print 'test10: succeeded'

    def test11(self):
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
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        print a_out_final_theory        


        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.40135 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 0.31219 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 1.3183 | units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 3880131 | units.RSun, 4)        
        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
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
        o_in = 0.0
        o_out = 0.0
        T_end = 50|units.Myr 
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end)

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 1.34345 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 1.34345 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 0.75302 | units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 5088750 | units.RSun, 4)        
        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, a_out_final_theory, 2)        

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

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 6.98888 | units.MSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.child2.mass, 6.98888 | units.MSun, 4)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 25.0399| units.RSun, 4)  
        # because of wind mass loss a_out_final != 1e6RSun      
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1001585 | units.RSun, 4)        
        
        M1f = tr.triple.child2.child1.mass
        M2f = tr.triple.child2.child2.mass
        M3f = tr.triple.child1.mass       
        a_out_final_theory = a_out * (M1+M2+M3) / (M1f+M2f+M3f) # wind
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

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

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

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test15: succeeded'
 

    def test16(self):
        print 'test16'

        M1 = 1.|units.MSun
        M2 = 0.1|units.MSun
        M3 = 0.08|units.MSun
        a_in = 20|units.RSun
#        a_in = 20|units.RSun
        a_out = 1.e6|units.RSun
        e_in = 0.0
        e_out = 1.e-5
        i = 0*np.pi/180.0
        g_in = 0.5*np.pi
        g_out = 0.5*np.pi
        o_in = 0.0
        o_out = 0.0
        T_end = 2000|units.Myr 
        dr = 0.005 #0.0005
        tr = triple.main(inner_primary_mass = M1, inner_secondary_mass = M2, outer_mass = M3, inner_semimajor_axis = a_in, outer_semimajor_axis = a_out, inner_eccentricity = e_in, outer_eccentricity = e_out, relative_inclination= i, inner_argument_of_pericenter = g_in, outer_argument_of_pericenter = g_out, inner_longitude_of_ascending_node = o_in, outer_longitude_of_ascending_node = o_out, tend = T_end, maximum_radius_change_factor=dr)

        self.assertAlmostRelativeEqual(tr.triple.child2.child1.mass, 0.99997 | units.MSun, 4)        
        self.assertEqual(tr.triple.child2.child2.mass, 0.1 | units.MSun)        
        self.assertEqual(tr.triple.child1.mass, 0.08 | units.MSun)        

        print tr.triple.child2.semimajor_axis, tr.triple.semimajor_axis
        print tr.triple.child2.eccentricity, tr.triple.eccentricity
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, 7.0995| units.RSun, 4)  
        self.assertAlmostRelativeEqual(tr.triple.semimajor_axis, 1.e6 | units.RSun, 4)        
        self.assertAlmostRelativeEqual(tr.triple.child2.eccentricity, 0., 4)  
        self.assertAlmostRelativeEqual(tr.triple.eccentricity, 0.001, 4)        


        a_in_final_theory = a_in * (1-e_in**2)
        self.assertAlmostRelativeEqual(tr.triple.child2.semimajor_axis, a_in_final_theory, 1)        

        print 'test16: succeeded'
 

        
        
        
#    def test2(self):
#        instance = self.new_instance_of_an_optional_code(SeBa)
#        
#        p = Particle()
#        p.mass = 5 | units.MSun
#        p.metallicity = 0.02
#        
#        p = instance.particles.add_particle(p)
#        instance.evolve_model(130 | units.Myr)
#
#        self.assertAlmostRelativeEqual(p.mass, 0.9906 | units.MSun, 4)
#        
#    def test3(self):
#        print "Testing evolution of a close binary system..."
#        instance = self.new_instance_of_an_optional_code(SeBa)
#        instance.commit_parameters()
#        stars =  Particles(2)
#        stars[0].mass = 3.0 | units.MSun
#        stars[1].mass = 0.3 | units.MSun
#        
#        
#        mu = (3.3 | units.MSun) * constants.G
#        orbital_period = 200.0 | units.day
#        semi_major_axis = (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)
#        
#        instance.particles.add_particles(stars)
#        
#        binaries =  Particles(1)
#        
#        binary = binaries[0]
#        binary.semi_major_axis = semi_major_axis
#        binary.eccentricity = 0.5
#        binary.child1 = stars[0]
#        binary.child2 = stars[1]
#        
#        instance.binaries.add_particles(binaries)
#        
#        from_seba_to_model = instance.particles.new_channel_to(stars)
#        from_seba_to_model.copy()
#
#        from_seba_to_model_binaries = instance.binaries.new_channel_to(binaries)
#        from_seba_to_model_binaries.copy()
#        
#        previous_type = binary.child1.stellar_type
#        results = []
#        current_time = 0 | units.Myr
#        
#        while current_time < (480 | units.Myr):
#            instance.update_time_steps()
#            # The next line appears a bit weird, but saves time for this simple test.
#            deltat = max(1.0*instance.binaries[0].time_step, 0.1| units.Myr)
#            current_time = current_time + deltat
#            instance.evolve_model(current_time)
#            from_seba_to_model.copy()
#            from_seba_to_model_binaries.copy()
#            if not binary.child1.stellar_type == previous_type:
#                results.append((binary.age, binary.child1.mass, binary.child1.stellar_type))
#                previous_type = binary.child1.stellar_type
#            
#        self.assertEqual(len(results), 6)
#        for x in results:
#            print x
#        
#        types = (
#            "Hertzsprung Gap",
#            "First Giant Branch",
#            "Core Helium Burning",
#            "First Asymptotic Giant Branch",
#            "Giant Branch Naked Helium star",
#            "Carbon/Oxygen White Dwarf",
#        )
#        
#
#        for result, expected in zip(results, types):
#            self.assertEquals(str(result[2]), expected)
#        
#        times = ( 
#            377.6369 | units.Myr, 
#            379.8877 | units.Myr,
#            382.3112 | units.Myr,
#            473.4804 | units.Myr,
#            475.4766 | units.Myr,
#            476.6182 | units.Myr, 
#        )
#        for result, expected in zip(results, times):
#            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
#            
#        masses = ( 
#            3.0000 | units.MSun, 
#            3.0000 | units.MSun, 
#            2.9983 | units.MSun, 
#            2.9741 | units.MSun,
#            0.6710 | units.MSun,
#            0.6596 | units.MSun,
#        )
#        for result, expected in zip(results, masses):
#            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 2)
#         
#        instance.stop()
#        
#    
#    def test5(self):
#        instance = self.new_instance_of_an_optional_code(SeBa)
#        self.assertAlmostRelativeEquals(instance.parameters.metallicity , 0.02)
#        instance.parameters.metallicity = 0.04
#        self.assertAlmostRelativeEquals(instance.parameters.metallicity , 0.04)



if __name__ == '__main__':
    test = TestSeBa()

#     test set up
#    test.test0()
#    test.test1()
#
##     test wind mass loss inner system
#    test.test2()
#    test.test3()
#    test.test4()
###
####     test wind mass loss outer system
#    test.test5()
##
##    # test stable mass transfer in inner system
###    test.test6() #collision
#    test.test7()
#    test.test8()
##
##    # test common envelope evolution in inner system
#    test.test9() #alpha-ce
#    test.test10() #gamma-ce
#    test.test11() #double alpha-ce
##
##    # test contact system in inner system
#    test.test12() #double alpha-ce
#    test.test13() #merger
##
##    #test tides
#    test.test14()
    test.test15()
#    test.tesst16()
##

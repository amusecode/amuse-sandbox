import pstats
p = pstats.Stats('prof')
#p.strip_dirs().sort_stats(-1).print_stats()
p.sort_stats('cumulative').print_stats()
#p.print_callees()

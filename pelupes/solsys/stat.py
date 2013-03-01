import pstats
p = pstats.Stats('prof')

p.strip_dirs().sort_stats('cumulative').print_stats(60)

p.strip_dirs().sort_stats('cumulative').print_callers()

p.strip_dirs().sort_stats('cumulative').print_callees()

import pstats
p = pstats.Stats('diomira.stat')
p.sort_stats('time').print_stats(10)
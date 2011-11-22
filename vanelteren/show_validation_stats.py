import pstats
import sys
import collections
import inspect

FunctionStatistics = collections.namedtuple(
    'FunctionStatistics', 
    'number_of_primitive_calls number_of_calls total_time_in_call cumulative_time callers'
)

def get_function_in_module(function_list, name_of_file, name_of_function, line_number):
    selected_functions = []
    for path, lineno, function_name in function_list:
        if path.endswith(name_of_file) and function_name == name_of_function:
            selected_functions.append( (path, lineno, function_name) )
    if len(selected_functions) == 1:
        return selected_functions[0]
        
    if len(selected_functions) == 0:
        raise Exception("could not find function {0} in file {1}".format(name_of_function, name_of_file))
    
    min_distance = 9999
    index_of_result = -1
    for i,x in enumerate(selected_functions):
        distance = abs(x[1] - line_number)
        if distance == 0:
            return x
        if distance < min_distance:
            index_of_result = i
            min_distance = distance
    print "warning requested line:", line_number, "but found function on line", selected_functions[index_of_result][1]
    return selected_functions[index_of_result]
    
    
    
def main(filename):
    stats = pstats.Stats(filename)
    stats.sort_stats('cumulative')
    stats.print_stats(20)
    
    stat_list = list(stats.fcn_list)
    
    module = stat_list[0]
    assert module[-1] == '<execfile>' or module[-1] == '<module>'
    format_per = "{0: >5.1f}%"
    format_sec = "{0: >7.1f}"
    module_stats = FunctionStatistics(*stats.stats[module])
    print "Total time:", format_sec.format(module_stats.cumulative_time)
    
    time_spent_in_codes = get_function_in_module(stat_list, 'amuse/rfi/channel.py', 'receive_header', 149)
    time_spent_in_codes_stats = FunctionStatistics(*stats.stats[evolve_model_in_bridge])
    print "In Code   :", format_sec.format(time_spent_in_codes_stats.cumulative_time), format_per.format(time_spent_in_codes_stats.cumulative_time/module_stats.cumulative_time * 100.0)
    receive_header
    if module[0] == '~' or module[0] == '<string>':
        nameformodule='particles_and_gas_in_cluster.py'
    else:
        nameformodule = module[0]
        
    
    evolve_model_func = get_function_in_module(stat_list, nameformodule, 'evolve_model', 289)
    evolve_model_stats = FunctionStatistics(*stats.stats[evolve_model_func])
   
    evolve_model_in_bridge = get_function_in_module(stat_list, 'amuse/couple/bridge.py', 'evolve_model', 459)
    evolve_model_stats = FunctionStatistics(*stats.stats[evolve_model_in_bridge])
    print "Evolve time :", format_sec.format(evolve_model_stats.cumulative_time)
    
    kick_codes_in_bridge = get_function_in_module(stat_list, 'amuse/couple/bridge.py', 'kick_codes', 459)
    kick_stats = FunctionStatistics(*stats.stats[kick_codes_in_bridge])
    print "Kick codes  :", format_sec.format(kick_stats.cumulative_time), format_per.format(kick_stats.cumulative_time / evolve_model_stats.cumulative_time * 100.0)
    
    drift_codes_in_bridge = get_function_in_module(stat_list, 'amuse/couple/bridge.py', 'drift_codes', 459)
    drift_stats = FunctionStatistics(*stats.stats[drift_codes_in_bridge])
    print "Drift codes :", format_sec.format(drift_stats.cumulative_time), format_per.format(drift_stats.cumulative_time / evolve_model_stats.cumulative_time * 100.0)
    
    #drift_codes_in_bridge = get_function_in_module(stat_list, 'amuse/couple/bridge.py', 'drift', 459)
    #drift_stats = FunctionStatistics(*stats.stats[drift_codes_in_bridge])
    #print "Drift in codes :", format_sec.format(drift_stats.cumulative_time), format_per.format(drift_stats.cumulative_time / evolve_model_stats.cumulative_time * 100.0)
    
    
    receive_overhead = get_function_in_module(stat_list, 'amuse/rfi/channel.py', 'receive_content', 459)
    receive_overhead_stats = FunctionStatistics(*stats.stats[receive_overhead])
    print "Receive     :",format_sec.format(receive_overhead_stats.cumulative_time)
    
    send_overhead = get_function_in_module(stat_list, 'amuse/rfi/channel.py', 'send_content', 459)
    send_overhead_stats = FunctionStatistics(*stats.stats[send_overhead])
    print "Send        :",format_sec.format(send_overhead_stats.cumulative_time)
    
    
    #stats.print_stats(20)
    #stats.print_stats('get_grav')
    #stats.print_stats('get_pot')
    #stats.print_stats('drift')
    #stats.print_stats('kick')
    # stats.print_stats('update_velocities')
    

if __name__ == '__main__':
    
    main(sys.argv[1])
    

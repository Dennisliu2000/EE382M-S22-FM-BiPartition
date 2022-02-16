# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 20:08:03 2022

@author: denni
"""

import os
from student_impl.eid_dl34437_EID import FM_Partition
eid = "dl34437"
benchmark_path = "benchmarks/example_2.txt"
output_root = "output"
output_root = os.path.join(output_root, eid)
if not os.path.isdir(output_root):
    os.mkdir(output_root)

output_path = os.path.join(output_root, os.path.basename(benchmark_path))
solver = FM_Partition()
solver.read_graph(benchmark_path)
solution = solver.solve()
profiling = (0, 0) # ignore runtime and memory for now
solver.dump_output_file(*solution, *profiling, output_path)
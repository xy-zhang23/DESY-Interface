/usr/gridengine/default/spool/bladek6/job_scripts/66033500:12: command not found: ^M
Traceback (most recent call last):
  File "injector_optimization_demo.py", line 37, in <module>
    pool = MPIPool()
  File "/afs/ifh.de/group/pitz/data/lixiangk/work/apps/python3/lib/python3.6/site-packages/platypus/mpipool.py", line 66, in __init__
    raise ValueError("Tried to create an MPI pool, but there "
ValueError: Tried to create an MPI pool, but there was only one MPI process available. Need at least two.

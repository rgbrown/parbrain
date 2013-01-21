# Example MPI LoadLeveler Job file
# @ shell = /bin/ksh
#
# @ job_name = parbrain
#
# @ job_type = parallel
#
# @ node = 8 
# @ tasks_per_node = 32  
#
# @ wall_clock_limit = 00:20:00
#
# Groups to select from: UC, UC_merit, NZ, NZ_merit
# @ group = UC_merit
#
# @ output = $(job_name).$(schedd_host).$(jobid).out
# @ error = $(job_name).$(schedd_host).$(jobid).err
# @ notification = never
# @ class = p7linux
#
# @ rset = rset_mcm_affinity
# @ task_affinity = core(1)
# @ network.MPI_LAPI = sn_single,shared,US,,instances=2
# @ queue
# suggested environment settings:
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MEMORY_AFFINITY=MCM

poe ./simulate 15 3

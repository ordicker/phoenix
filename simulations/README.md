# Installing MEEP 

- install conda/miniconda
- make venv \w parallel meep
```console
conda create -n meep -c conda-forge pymeep=*=mpi_mpich_*
```

# Activating meep venv 
```console 
conda activate meep
```

# Installing pylsp
```console
conda install -c conda-forge python-lsp-server
```

# running script with 4 threads
```console
mpirun -np 4 python <script_name>.py
```

# problem solving 

## gsl version error 
got this error 
```console
ImportError: libgsl.so.25: cannot open shared object file: No such file or directory                                              
```
reverting from gsl=2.7.1 to gsl=2.7.0 solved the problem
```console
conda install gsl=2.7.0
```


gfortran  -g -fbacktrace -Wall -fcheck=all -cpp -DUSE_NETCDF\
    src/teco/datatypes.f90 \
    src/teco/io_mod.f90 \
    src/teco/soil.f90 \
    src/teco/vegetation.f90 \
    src/teco/transfer.f90 \
    src/tools/mcmc_mod.f90 \
    src/tools/mcmc_outputs.f90  \
    src/teco/driver.f90 \
    src/tools/mcmc.f90   \
    `# src/tools/spinup.f90`  \
    src/main.f90 \
    -o run_teco \
    -I/home/zhou_j/miniconda3/envs/spruce-da/include -L/home/zhou_j/miniconda3/envs/spruce-da/lib -lnetcdff -lnetcdf

if find "src" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/*.mod
fi

if find "src/teco" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/teco/*.mod
fi

if find "src/tools" -name "*.mod" -print -quit | grep -q '.*'; then
    rm src/tools/*.mod
fi

if find "." -name "*.mod" -print -quit | grep -q '.*'; then
    rm *.mod
fi

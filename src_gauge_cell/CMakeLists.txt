cmake_minimum_required(VERSION 3.16)
project(src_gauge_cell C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(src_gauge_cell
cubic_spline_1d.h
        framework_force.c    integration.h       inter_energy.c      mc_moves.c                 molecule_properties.h  recrossing.c          simulation.h         utils.c
         equations_of_state.c  framework_force.h    internal_born.c     inter_energy.h      mc_moves.h                 monte_carlo.c          recrossing.h          spacegroup.c         utils.h
cbmc.c                  equations_of_state.h  framework_hessian.c  internal_born.h     inter_force.c       minimization.c             monte_carlo.h          rigid.c               spacegroup.h         vector.c
cbmc.h                  ewald.c               framework_hessian.h  internal_energy.c   inter_force.h       minimization.h             movies.c               rigid.h               spectra.c            vector.h
charge_equilibration.c  ewald.h               grids.c              internal_energy.h   inter_hessian.c     minimization_simulation.c  movies.h               run.c                 spectra.h            warnings.c
charge_equilibration.h  framework.c           grids.h              internal_force.c    inter_hessian.h     minimization_simulation.h  numerical.c            run.h                 statistics.c         warnings.h
complex.c               framework.h           input.c              internal_force.h    linear_equations.c  molecular_dynamics.c       numerical.h            sample.c              statistics.h
complex.h               framework_born.c      input.h              internal_hessian.c  linear_equations.h  molecular_dynamics.h       output.c               sample.h              status.c
constants.c             framework_born.h      integrate.c          internal_hessian.h  main.c              molecule.c                 output.h               scattering_factors.c  status.h
constants.h             framework_energy.c    integrate.h          inter_born.c        matrix.c            molecule.h                 potentials.c           scattering_factors.h  thermo_baro_stats.c
cubic_spline_1d.c       framework_energy.h    integration.c        inter_born.h        matrix.h            molecule_properties.c      potentials.h           simulation.c          thermo_baro_stats.h
        gauge_cell.c gauge_cell.h)

target_link_libraries(src_gauge_cell m)

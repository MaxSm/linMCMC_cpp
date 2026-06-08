cmake_minimum_required(VERSION 3.30)

project(solver)

add_executable(lin_ode_mcmc src/lin_ode_mcmc.cpp src/lin_ode_mcmc_functions.cpp src/lin_ode_mcmc.h)
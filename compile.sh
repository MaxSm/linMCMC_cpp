#!/bin/bash
g++ -c lin_ode_mcmc.cpp
g++ -c lin_ode_mcmc_functions.cpp
g++ -o lin_ode_mcmc lin_ode_mcmc.o lin_ode_mcmc_functions.o

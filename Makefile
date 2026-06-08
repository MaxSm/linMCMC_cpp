# Targets
solver: lin_ode_mcmc

solver_db: lin_ode_mcmc_db

# Link object files 
lin_ode_mcmc: lin_ode_mcmc.cpp lin_ode_mcmc_functions.cpp
	g++ -o lin_ode_mcmc lin_ode_mcmc.cpp lin_ode_mcmc_functions.cpp

lin_ode_mcmc_db: lin_ode_mcmc.cpp lin_ode_mcmc_functions.cpp
	g++ -g -o lin_ode_mcmc_db lin_ode_mcmc.cpp lin_ode_mcmc_functions.cpp


# Clean up built files
clean: 
	rm lin_ode_mcmc
clean_db:
	rm lin_ode_mcmc_db
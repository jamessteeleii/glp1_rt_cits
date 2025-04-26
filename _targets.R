# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed.
library(crew)

# Set target options:
tar_option_set(packages = c("here",
                            "tidyverse",
                            "metafor",
                            "faux",
                            "nlme",
                            "marginaleffects",
                            "furrr",
                            "progressr"),
               seed = 1988,  # <-- GLOBAL reproducible seed
               memory = "transient",
               format = "qs",
               garbage_collection = TRUE,
               storage = "worker",
               retrieval = "worker",
               controller = crew_controller_local(workers = 20)) # Packages that your targets need for their tasks.
               # format = "qs", # Optionally set the default storage format. qs is fast.
               #
               # Pipelines that take a long time to run may benefit from
               # optional distributed computing. To use this capability
               # in tar_make(), supply a {crew} controller
               # as discussed at https://books.ropensci.org/targets/crew.html.
               # Choose a controller that suits your needs. For example, the following
               # sets a controller that scales up to a maximum of two workers
               # which run as local R processes. Each worker launches when there is work
               # to do and exits if 60 seconds pass with no tasks to run.
               #
               #   controller = crew::crew_controller_local(workers = 2, seconds_idle = 60)
               #
               # Alternatively, if you want workers to run on a high-performance computing
               # cluster, select a controller from the {crew.cluster} package.
               # For the cloud, see plugin packages like {crew.aws.batch}.
               # The following example is a controller for Sun Grid Engine (SGE).
               #
               #   controller = crew.cluster::crew_controller_sge(
               #     # Number of workers that the pipeline can scale up to:
               #     workers = 10,
               #     # It is recommended to set an idle time so workers can shut themselves
               #     # down if they are not running tasks.
               #     seconds_idle = 120,
               #     # Many clusters install R as an environment module, and you can load it
               #     # with the script_lines argument. To select a specific verison of R,
               #     # you may need to include a version string, e.g. "module load R/4.3.2".
               #     # Check with your system administrator if you are unsure.
               #     script_lines = "module load R"
               #   )
               #
               # Set other options as needed.)
               
               # Run the R scripts in the R/ folder with your custom functions:
               tar_source("R/functions.R")
               # tar_source("other_functions.R") # Source other scripts as needed.
               
               # Replace the target list below with your own:
               list(
                 ##### Pre-registration targets ----
                 
                 ### Read in and prepare data to use for determining parameters in simulations
                 
                 # Load and prepare data from previous meta-analyses and studies
                 tar_target(karakakis_data_file,
                            here("data", "karakakis_baseline_study_data.csv"),
                            format = "file"),
                 tar_target(karakakis_data,
                            read_prep_karakakis_data(karakakis_data_file)),
                 
                 tar_target(benito_data_file,
                            here("data", "benito_study_data.csv"),
                            format = "file"),
                 tar_target(benito_data,
                            read_prep_benito_data(benito_data_file)),
                 
                 tar_target(murphy_koehler_data_file,
                            here("data", "murphy_koehler_study_data.csv"),
                            format = "file"),
                 tar_target(murphy_koehler_data,
                            read_prep_murphy_koehler_data(murphy_koehler_data_file)),
                 
                 tar_target(krupa_data, prep_krupa_data()),
                 
                 tar_target(beaver_data_file,
                            here("data", "beaver_study_data.csv"),
                            format = "file"),
                 tar_target(beaver_data, read_prep_beaver_data(beaver_data_file)),
                 
                 ### Prepare parameters for simulations and check assumptions
                 tar_target(raw_glp1_effect, get_raw_glp1_effect()),
                 
                 tar_target(SMD_glp1_effect,
                            get_SMD_glp1_effect(karakakis_data, raw_glp1_effect)),
                 
                 tar_target(raw_RT_effect,
                            get_raw_RT_effect(benito_data)),
                 
                 tar_target(SMD_RT_effect,
                            get_SMD_RT_effect(benito_data)),
                 
                 tar_target(SMD_RT_energy_restriction_effect,
                            get_SMD_RT_energy_restriction_effect(murphy_koehler_data)),
                 
                 tar_target(SMD_energy_restriction_effect,
                            get_SMD_energy_restriction_effect(krupa_data)),
                 
                 # Check assumption of interindividual response variation to GLP1-RA on lean mass
                 tar_target(irv_glp1_effect,
                            check_glp1_irv_assumption(beaver_data)),
                 
                 tar_target(intercepts_lean_mass,
                            set_intercepts_lean_mass()),
                 tar_target(measurement_error,
                            set_measurement_error()),
                 
                 
                 ### Set parameters and run simulation for power
                 tar_target(simulation_parameters,
                            define_simulation_parameters(intercepts_lean_mass,
                                                         raw_RT_effect,
                                                         raw_glp1_effect,
                                                         measurement_error)),
                 
                 tar_target(
                   simulation_parameter_batches,
                   split(simulation_parameters, ceiling(row_number(simulation_parameters)/1000)),
                   iteration = "list"  # important to pass list elements separately
                 ),
                 
                 tar_target(
                   simulation_results,
                   simulate_power(simulation_parameter_batches),
                   pattern = map(simulation_parameter_batches)
                  ),
                 
                 tar_target(
                   simulation_results_combined,
                   bind_rows(simulation_results)
                 )
                 
                 # tar_target(power_simulations,
                 #            simulate_power(simulation_parameters))

                 
               )
               
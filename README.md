## Overview

This R project contains the code required to replicate the simulations in the current revision of ENHANCE (E-value Nested & Hierarchical Analysis of Networks, Controlling Error). In order to run the simulations, the abind and CVXR libraries are required (and therefore, so is an up-to-date Rust installation), while to run the visualizations, the ggplot and gt libraries are required.

The simulations take between 3 and 6 hours to run in full and are memory-intensive. To aid in preventing loss of work, the simulation process is split into three steps: the generation of the simulated data, which takes a short time; the testing procedure, which is the most time-intensive segment; and the visualization step; across two files, one for simulations and one for visualizations. Outputs are saved after the first two steps to prevent loss of work. **You may be required to make your own outputs/ subdirectory before running any simulations. This is where saved outputs are directed by default.**

## Files

e_procedures.R: bearing a slight misnomer, this file contains the p- and e-value elicitation procedures.

sim_utils.R: this file contains the functions required to carry out the simulation step.

viz_utils.R: this file contains the functions required to carry out the visualization step.

sims.R: this file contains code to run the simulations.

viz.R: this file contains code to run the visualizations.

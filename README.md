# Phase-locked TMS for Essential Tremor Simulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the MATLAB/NEURON simulation code for our research on phase-locked transcranial magnetic stimulation (TMS) for essential tremor treatment.

## Publication

Zhang, X., Hancock, R., Santaniello, S. (2025) Feasibility of Phase-Locked Transcranial Magnetic Stimulation of Cerebellum for the Treatment of Essential Tremor, *Journal of Neural Engineering*.

## Prerequisites

- **MATLAB** (recommended: R2020a or newer)
- **NEURON** (recommended: 8.0.0 or newer)

## Repository Structure

The repository contains two main simulation directories:
- `PRC_scaleup_5x_6d3Hz`: Simulations for 6.3 Hz tremor frequency
- `PRC_scaleup_5x_7d2Hz`: Simulations for 7.2 Hz tremor frequency

## Setup Instructions

1. **Compile MOD files**:
   - Navigate to either tremor frequency directory
   - Compile the MOD files in the `modfiles` folder using `mknrndll`
   - This will generate `nrnmech.dll` (Windows) or an `x86_64` folder (Unix)
   - Move the generated files to the parent directory (outside the `modfiles` folder)

2. **Initialize Simulation**:
   - Run `ET_start.hoc` to initialize the network states (simulates first 2,000 ms)
   - This generates state files in the `initialization` folder needed for subsequent simulations

## Running Simulations

After initialization, you can run any of the following MATLAB scripts:

### TMS Simulations

| Script | Description |
|--------|-------------|
| `run_rTMS_par_*p.m` | Regular/periodic rTMS (open-loop TMS) at various frequencies |
| `run_irTMS_par_*p.m` | Irregularly spaced TMS pulses (Sobol sequence) |
| `run_TBS_par_*p.m` | Continuous theta-burst TMS |
| `run_PL_TMS_par_*p_s*.m` | Phase-locked TMS with s* cycles skipped between pulses |

*Note: `*p` indicates the percentage of Purkinje cells activated by each TMS pulse*

#### Reproducing Figure 5

To reproduce Figure 5 from the paper:

1. Navigate to the `PRC_scaleup_5x_7d2Hz` directory
2. Run `run_PL_TMS_par_20p_s0.m`
3. For simplicity, only two simulations need to be run (one for an effective phase and another for an ineffective phase)
4. Modify the script by changing:
   ```matlab
   run_range = 1:length(PRClist);
   ```
   to:
   ```matlab
   run_range = [12 46];
   ```
5. Run `Figure_5.m` to obtain the plots

### tACS Simulations

| Script | Description |
|--------|-------------|
| `run_OL_tACS.m` | Open-loop tACS |
| `run_PL_tACS_*pA.m` | Phase-locked tACS (`*pA` indicates current magnitude) |

All simulations run until 10,000 ms.

## Optional Configuration

To create differently randomized synaptic connections, run `createConnections.m` in MATLAB with a different `rngSEED` value.

## Contact

For questions regarding simulations or analyses, please contact:  
Xu Zhang: Xu.Zhang@childrens.harvard.edu
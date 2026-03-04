# Windkessel Parameter Estimation (MATLAB)

MATLAB code for estimating **three-element Windkessel (WK3)** parameters from an input flow waveform by solving the governing ODE and running an optimization.

For full documentation, background, and theory, please see the project page:
https://zongze-li-pon.github.io/completed-projects/2023-06-15-windkessel-parameter-estimation-matlab-code/

---

## Quick Start

1. **Prepare input data**  
   Put your waveform into `flowrate.csv` (Column 1 = time, Column 2 = flow rate).

2. **Run the script in MATLAB**  
   Open MATLAB in this folder and run:
  
   windkessel_model_parameters_estimation

3. **Check outputs**  
The script prints the estimated parameters and generates the pressure waveform plot.

---

## More Details

The full methodology and background of this work are described in the paper:

A fast approach to estimating Windkessel model parameters for patient-specific multi-scale CFD simulations of aortic flow  
https://www.sciencedirect.com/science/article/abs/pii/S0045793023001196
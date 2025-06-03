# Robust Adaptive Model Predictive Control for Tracking in Interconnected Systems via Distributed Optimisation

This repository contains all the files required to test the DRAMPC algorithm presented in the paper for the 5 agents setup.

*Note that the scripts can be modified to simulate up to 9 agents, as the toolbox MPT3 cannot manage higher dimension polyhedra in the centralised setup*

## Repository Structure

- **Classes/**  
  Contains all the necessary MATLAB classes to run the code. 
  - `MPC`: Model Predictive Control
  - `RMPC`: Robust MPC -> Tube-based MPC
  - `RAMPC`: Robust Adaptive MPC
  - `DMPC, DRMPC, DRAMPC`: Distributed MPC flavours.
  
  The class hierarchy is as follows:  
  - `MPC < RMPC < RAMPC`  
  - `MPC < DMPC < DRMPC < DRAMPC`  
  where `<` denotes class inheritance.  
  Additional classes include:
  - `Agents`: Defines the agents and their linear dynamics.
  - `Network`: Defines the network topology in which the agents operate.

- **Scripts/**  
  Contains scripts to test the algorithms. Three scripts are provided for five agents, each testing MPC, DMPC, RAMPC, and DRAMPC.  
  The scripts test:
  - Tracking of piecewise constant setpoints
  - Tracking of a piecewise constant signal trajectory
  - Tracking of a harmonic signal trajectory

- **Utils/**  
  Contains utility functions used by the scripts and classes.

- **Main Directory**  
  - `requirements.txt`: List of required MATLAB toolboxes and packages.
  - `install_req.m`: MATLAB script to install missing requirements. (BETA available soon)

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](./LICENSE) file for details.



# BeSim Toolbox
Matlab toolbox for quick design and simulation of advanced building climate control algorithms.

## Features
- Interface for [linearized white-box building envelope models from Modelica](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- Automated construction and tuning of model predictive control (MPC) and state estimation algorithms
- Closed-loop simulation, plotting, and performance analysis
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903) (deep lerarning in particular)
- FMI interface for Modelica emulator models (comming soon)
- For quick start and more details check out the presentation about algorithms and tools behind [BuiSim](https://www.researchgate.net/publication/328171184_Tools_and_Techniques_for_Advanced_Model_Predictive_Building_Control)

## Installation
### tbxmanager
- Install [tbxmanager](http://www.tbxmanager.com/) 
- Install BeSim via:  `tbxmanager install besim`
- Check for updates:  `tbxmanager update` 
### manual
- clone  BeSim repository 
- save BeSim folder with its subfolders to Matlab path  

## Prerequisities
- Matlab: developed and tested on R2017a and R2017b
- [Yalmip](https://yalmip.github.io/) mathematical modeling and optimization toolbox (BeSim's backbone)
- Advanced optimization solver, e.g. [Gurobi](http://www.gurobi.com/) (solution of implicit MPC and MHE problems)
- Matlab toolboxes: Deep Learning, Machine learning (approximate MPC functionality)

## Getting Started - Demos
run following scripts in Matlab to get quick results:
- [BeInit.m](https://github.com/drgona/BeSim/blob/master/Be_Run/BeInit.m): design and simulation of optimization-based MPC and state estimator for selected building model
- [BeInitML.m](https://github.com/drgona/BeSim/blob/master/Be_Run/BeInitML.m): design and simulation of approximate MPC via machine learning for selected buiding model

## Structure
**Functional Structure:** Graphical overview of BuiSim structure with data-flow dependencies.
![BuiSim structure](/Data/Page/BeSim_structure2.png)

**Repository Structure:**
List of repository folders with associated functionality.
- [Be_Run](https://github.com/drgona/BeSim/tree/master/Be_Run) (run files and demos)
- [Be_Modeling](https://github.com/drgona/BeSim/tree/master/Be_Modeling) (model loading, discretization, model order reduction)
- [Be_Disturbances](https://github.com/drgona/BeSim/tree/master/Be_Disturbances) (disturbance trajectories loading)
- [Be_References](https://github.com/drgona/BeSim/tree/master/Be_References) (reference trajectories loading)
- [Be_Estimation](https://github.com/drgona/BeSim/tree/master/Be_Estimation) (estimator design)
- [Be_Control](https://github.com/drgona/BeSim/tree/master/Be_Control) (controller design)
- [Be_Simulation](https://github.com/drgona/BeSim/tree/master/Be_Simulation) (main simulation and plotting functions)
- [Be_Learn](https://github.com/drgona/BeSim/tree/master/Be_Learn) (machine learning functions, synthesis of approximate MPC)
- [buildings](https://github.com/drgona/BeSim/tree/master/buildings) (building models files)
- [Data](https://github.com/drgona/BuiSim/tree/master/Data) (stored results)

## Algorithms 
List of key enabling algorithms implemented in BeSim.

**Model Order Reduction**
- [Balanced truncation](https://nl.mathworks.com/help/robust/ref/reduce.html)

**State Estimation**
- Kalman filters (KF)
- Moving horizon estimation (MHE)

**Optimal Control**
- Model predictive control (MPC)
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903)

**Machine Learning Models**
- Deep learning (DL)
- Regression trees (RT)

## Building Models

**Model Structure**
- [Linearized white-box building envelope models](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- arbitrary linear time invariant model in state space format

**Available Building Models**

Building type | Location      |  Label        | floor area [m2] | #states         | #outputs       | #inputs         | #disturbances
------------  | ------------- | ------------- | -------------   | -------------  | -------------   | -------------  | ------------- 
Residential   |  Belgium      | 'Old', 'Reno', 'RenoLight' | 56 | 283,286,250    | 6               | 6               | 44
Office        |  Belgium      | 'HollandschHuys' | 3760         | 700            | 12              | 73               | 289
Office        |  Belgium      | 'Infrax'         | 2232         | 1262           | 19              | 28               | 259
Borehole      |  Belgium      | 'Borehole '      | -            | 190            | 1               | 1                   | 0

## Contact
Email: jan.drgona@kuleuven.be 

Author: [Ján Drgoňa](https://www.kuleuven.be/wieiswie/en/person/00107194)  
postdoctoral researcher  
KU Leuven  
Department of Mechanical Engineering  
Division of Applied Mechanics and Energy Conversion (TME)  
Celestijnenlaan 300A, BE-3001 Leuven (Heverlee), Belgium  





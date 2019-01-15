# BuiSim Toolbox
Matlab toolbox for quick design and simulation of advanced building climate control algorithms.

## Features
- Interface for [linearized white-box building envelope models from Modelica](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- Automated construction and tuning of model predictive control (MPC) and state estimation algorithms
- Closed-loop simulation, plotting, and performance analysis
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903) (deep lerarning in particular)
- FMI interface for Modelica emulator models (comming soon)
- For quick start and more details check out the presentation about algorithms and tools behind [BuiSim](https://www.researchgate.net/publication/328171184_Tools_and_Techniques_for_Advanced_Model_Predictive_Building_Control)

## Installation
- clone the BuiSim repository
- save BuiSim folder with its subfolders to Matlab path  

## Prerequisities
- [Yalmip](https://yalmip.github.io/) mathematical modeling and optimization toolbox: BuiSim's backbone.
- Advanced optimization solver, e.g. [Gurobi](http://www.gurobi.com/): for solution of implicit MPC, MHE problems.
- Matlab toolboxes: Deep Learning, Machine learning (approximate MPC functionality)

## Getting Started - Demos
run following scripts in Matlab to get quick results:
- [BuiInit.m](https://github.com/drgona/BuiSim/blob/master/Bui_Run/BuiInit.m): design and simulation of optimization-based MPC and state estimator for selected building model
- [BuiInitML.m](https://github.com/drgona/BuiSim/blob/master/Bui_Run/BuiInitML.m): design and simulation of approximate MPC via machine learning for selected buiding model

## Structure
**Functional Structure:** Graphical overview of BuiSim structure with data-flow dependencies.
![BuiSim structure](/Data/Page/BuiSim_structure2.png)

**Repository Structure:**
List of repository folders with associated functionality.
- [Bui_Run](https://github.com/drgona/BuiSim/tree/master/Bui_Run) (run files and demos)
- [Bui_Modeling](https://github.com/drgona/BuiSim/tree/master/Bui_Modeling) (model loading, discretization, model order reduction)
- [Bui_Disturbances](https://github.com/drgona/BuiSim/tree/master/Bui_Disturbances) (disturbance trajectories loading)
- [Bui_References](https://github.com/drgona/BuiSim/tree/master/Bui_References) (reference trajectories loading)
- [Bui_Estimation](https://github.com/drgona/BuiSim/tree/master/Bui_Estimation) (estimator design)
- [Bui_Control](https://github.com/drgona/BuiSim/tree/master/Bui_Control) (controller design)
- [Bui_Simulation](https://github.com/drgona/BuiSim/tree/master/Bui_Simulation) (main simulation and plotting functions)
- [Bui_Learn](https://github.com/drgona/BuiSim/tree/master/Bui_Learn) (machine learning functions, synthesis of approximate MPC)
- [buildings](https://github.com/drgona/BuiSim/tree/master/buildings) (building models files)
- [Data](https://github.com/drgona/BuiSim/tree/master/Data) (stored results)

## Algorithms 
List of key enabling algorithms implemented in BuiSim.

**Optimal Control**
- Model predictive control (MPC)
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903)

**State Estimation**
- Kalman filters
- Moving horizon estimation

**Machine Learning Models**
- Deep learning
- Regression trees

## Building Models

**Model Structure**
- [Linearized white-box building envelope models](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- arbitrary linear time invariant model in state space format

**Available Building Models**

Building type | Location      |  Label        | floor area [m2] | #states         | #outputs       | #inputs         | #disturbances
------------  | ------------- | ------------- | -------------   | -------------  | -------------   | -------------  | ------------- 
Residential   |  Belgium      | 'Old', 'Reno', 'RenoLight'  | 56 | 283,286,250 | 6    | 6               | 44
Office   |  Hasselt, Belgium      | 'HollandschHuys' | 3760 | 700 | 12    | 73               | 289
Office   |  Belgium      | 'Infrax' | 2232 | 1262 | 19    | 28               | 259

## Contact
- Send feedback to: jan.drgona@kuleuven.be 
- Author: [Ján Drgoňa](https://www.kuleuven.be/wieiswie/en/person/00107194)
postdoctoral researcher
Department of Mechanical Engineering
Division of Applied Mechanics and Energy Conversion (TME)
KU Leuven
Celestijnenlaan 300A, BE-3001 Leuven (Heverlee), Belgium





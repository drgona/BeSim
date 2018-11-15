
For quick start and more details check out the presentation about algorithms and tools behind [BuiSim](https://www.researchgate.net/publication/328171184_Tools_and_Techniques_for_Advanced_Model_Predictive_Building_Control).


## Features
- Matlab toolbox for quick design and simulation of advanced building climate control algorithms
- Interface for [linearized white-box building envelope models from Modelica](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- Automated construction and tuning of model predictive control (MPC) and state estimation algorithms
- Closed-loop simulation, plotting, and performance analysis
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903) (deep lerarning in particular)
- FMI interface for Modelica emulator models (comming soon)


## Structure
**Functional Structure**
![BuiSim structure](/Data/Page/BuiSim_structure2.png)

**Repository Structure:**
List of repository folders with associated functionality.
- Bui_Run (run files)
- Bui_Modeling (model loading, discretization, model order reduction)
- Bui_Disturbances (disturbance trajectories loading)
- Bui_References (reference trajectories loading)
- Bui_Estimation (estimator design)
- Bui_Control (controller design)
- Bui_Simulation (main simulation and plotting functions)
- Bui_Learn (machine learning functions, synthesis of approximate MPC)
- buildings (building models files)
- Data (stored results)

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


## Prerequisities
- [Yalmip](https://yalmip.github.io/) mathematical modeling and optimization toolbox: BuiSim's backbone.
- Advanced optimization solver, e.g. [Gurobi](http://www.gurobi.com/): for solution of implicit MPC, MHE problems.
- Matlab toolboxes, mainly Deep Learning, Machine learning: for approximate MPC functionality





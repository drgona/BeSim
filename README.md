
## Algorithms and Models
List of key enabling algorithms and models implemented in BuiSim.

**Building Models**
- [Linearized white-box building envelope models](http://www.ep.liu.se/ecp/article.asp?issue=118&article=005&volume=)
- arbitrary linear time invariant model in state space format

**Optimal Control**
- Model predictive control (MPC)
- [Approximate MPC via machine learning](https://www.sciencedirect.com/science/article/pii/S0306261918302903)

**State Estimation**
- Kalman filters
- Moving horizon estimation

**Machine Learning Models**
- Deep learning
- Regression trees

## Features
- Matlab interface for loading linearized white-box building envelope models from Modelica
- Automated construction and tuning of MPC and state estimation algorithms
- Closed-loop simulation, plotting, and performance analysis
- Approximate MPC via machine learning (deep lerarning in particular)
- FMI interface for Modelica emulator models (comming soon)

## Structure
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

## Prerequisities
- [Yalmip](https://yalmip.github.io/)
- Advanced optimization solver, e.g. [Gurobi](http://www.gurobi.com/)
- Matlab toolboxes: Deep Learning, Machine learning

```
![Image](src)

## Welcome to GitHub Pages

You can use the [editor on GitHub](https://github.com/drgona/BuiSim/edit/master/README.md) to maintain and preview the content for your website in Markdown files.

Whenever you commit to this repository, GitHub Pages will run [Jekyll](https://jekyllrb.com/) to rebuild the pages in your site, from the content in your Markdown files.


For more details see [GitHub Flavored Markdown](https://guides.github.com/features/mastering-markdown/).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/drgona/BuiSim/settings). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://help.github.com/categories/github-pages-basics/) or [contact support](https://github.com/contact) and we’ll help you sort it out.

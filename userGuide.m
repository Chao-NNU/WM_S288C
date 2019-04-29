% User guide

%% Lesson 1: generate the knowledgeBase data from MySQL Database
%go into the root directory of whole-cell toolbox
setWarnings();
setPath();
install%follow the instructions on screen
[knowledgeBase, knowledgeBaseWID] = cacheKnowledgeBase();

%% lesson 2: construct simulation
% import classes
import edu.jiangnan.fmme.cell.sim.Simulation;
import edu.jiangnan.fmme.cell.sim.util.FitConstants;
import edu.jiangnan.fmme.cell.sim.util.CachedSimulationObjectUtil;

%% recalculate knowledge base properties
knowledgeBase.invalidate();
knowledgeBase.calcIndices();
knowledgeBase.computeDependentProperties();

% generate simulation
simulation = Simulation(knowledgeBase.states, knowledgeBase.processes);

%initialize its constants
simulation.compartment.initializeConstants(knowledgeBase,simulation);
simulation.gene.initializeConstants(knowledgeBase,simulation);

% write simulation data
save('data/Simulation_fitted.mat', 'simulation', 'knowledgeBaseWID');

%initialize constants for states and processes
simulation.initializeConstants(knowledgeBase)
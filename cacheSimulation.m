% Caches simulation object.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/26/2012
function simulation = cacheSimulation(knowledgeBase, knowledgeBaseWID)
import edu.jiangnan.fmme.cell.sim.Simulation;
import edu.jiangnan.fmme.cell.sim.util.FitConstants;
import edu.jiangnan.fmme.cell.sim.util.CachedSimulationObjectUtil;

%% recalculate knowledge base properties
knowledgeBase.invalidate();
knowledgeBase.calcIndices();
knowledgeBase.computeDependentProperties();

%% construct simulation and initialize its constants
simulation = Simulation(knowledgeBase.states, knowledgeBase.processes);
simulation.initializeConstants(knowledgeBase);
simulation.applyOptions('seed', 1);
fitter = FitConstants(simulation, struct('verbosity', 1));
fitter.run();

%initializeConstants(simulation.compartment, knowledgeBase);
%initializeConstants(simulation.gene, knowledgeBase);

%% choose an initial simulation state with desired growth rate
simulation.applyOptions('seed', 1);
mr = simulation.state('MetabolicReaction');
initialGrowthFilterWidth = mr.initialGrowthFilterWidth;
mr.initialGrowthFilterWidth = 1e-2;
simulation.initializeState();
mr.initialGrowthFilterWidth = initialGrowthFilterWidth;

%% cache simulation object
CachedSimulationObjectUtil.store(simulation, knowledgeBaseWID);
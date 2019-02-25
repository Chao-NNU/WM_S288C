% Runs the medium tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runMediumTests()
%% initialize
warning('off','WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.jiangnan.fmme.test.XMLTestRunDisplay;
import edu.jiangnan.fmme.test.runtests;
monitor = XMLTestRunDisplay('Medium Tests','Whole cell simulation medium tests', 'output/runMediumTests/tests.xml');
runtests(monitor, {
    'edu.jiangnan.fmme.cell.sim.CellCycle_Test'
    'edu.jiangnan.fmme.cell.sim.Cytokinesis_Test'
    'edu.jiangnan.fmme.cell.sim.DNADamageRepair_Test'
    'edu.jiangnan.fmme.cell.sim.ProteinGrowth_Test'
    'edu.jiangnan.fmme.cell.sim.RNA_Test'
    'edu.jiangnan.fmme.cell.sim.Simulation_Integrated_Test'
    'edu.jiangnan.fmme.cell.sim.Simulation_KnowledgeBase_Test'
    'edu.jiangnan.fmme.cell.sim.Simulation_Test'
    'edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect_Test'    
    });


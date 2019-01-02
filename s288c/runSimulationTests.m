% Runs a full length simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/12/2011
function runSimulationTests(runProfiler)
%% profile on
if nargin >= 1 && runProfiler
    profile('on', '-nohistory');
end

%% initialize
setWarnings();
setPath();
setPreferences();

%% run tests
monitor = edu.jiangnan.fmme.test.XMLTestRunDisplay('Simulation Test','Whole cell full length simulation test', 'output/runSimulationTests/tests.xml');
edu.jiangnan.fmme.test.runtests(monitor, {
    'edu.jiangnan.fmme.cell.sim.Simulation_FullLength_Test'
    });

%% profile off, save data
if nargin >= 1 && runProfiler
    profile('off');
    profData = profile('info'); %#ok<NASGU>
    save('output/runSimulationTests/profData.mat', 'profData');
    edu.jiangnan.fmme.cell.sim.analysis.RunTime.save('output/runSimulationTests/profSummary.html');
end

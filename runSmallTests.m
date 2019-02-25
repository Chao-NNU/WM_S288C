% Runs the small tests for the whole cell simulation.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010
function runSmallTests()
%% initialize
warning('off', 'WholeCell:warning');
setWarnings();
setPath();
setPreferences();

%% run tests
import edu.jiangnan.fmme.test.XMLTestRunDisplay;
import edu.jiangnan.fmme.test.runtests;
monitor = XMLTestRunDisplay('Small Tests', 'Whole cell simulation small tests', 'output/runSmallTests/tests.xml');
runtests(monitor, {
    'edu.jiangnan.fmme.cell.sim.constant'
    'edu.jiangnan.fmme.cell.sim.process'
    'edu.jiangnan.fmme.cell.sim.state'
    'edu.jiangnan.fmme.cell.sim.util'
    'edu.jiangnan.fmme.cell.kb'
    'edu.jiangnan.fmme.db'
    'edu.jiangnan.fmme.io'
    'edu.jiangnan.fmme.test'
    'edu.jiangnan.fmme.util'
    });

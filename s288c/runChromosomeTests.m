% Run chromosome tests.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/13/2010

%% initialize
setWarnings();
setPath();
setPreferences();

%% run tests
runtests('-verbose', {
    'edu.jiangnan.fmme.cell.sim.state.Chromosome_Test'
    'edu.jiangnan.fmme.cell.sim.process.ChromosomeCondensation_Test'
    'edu.jiangnan.fmme.cell.sim.process.ChromosomeSegregation_Test'
    'edu.jiangnan.fmme.cell.sim.process.Cytokinesis_Test'
    'edu.jiangnan.fmme.cell.sim.process.DNADamage_Test'
    'edu.jiangnan.fmme.cell.sim.process.DNARepair_Test'
    'edu.jiangnan.fmme.cell.sim.process.DNASupercoiling_Test'
    'edu.jiangnan.fmme.cell.sim.process.Replication_Test'
    'edu.jiangnan.fmme.cell.sim.process.ReplicationInitiation_Test'
    'edu.jiangnan.fmme.cell.sim.process.Transcription_Test'
    'edu.jiangnan.fmme.cell.sim.process.TranscriptionalRegulation_Test'
    'edu.jiangnan.fmme.cell.sim.DNA_Test'
    'edu.jiangnan.fmme.cell.sim.DNADamageRepair_Test'
    });
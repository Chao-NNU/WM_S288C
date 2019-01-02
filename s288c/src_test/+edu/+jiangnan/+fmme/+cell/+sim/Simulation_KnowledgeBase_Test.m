%Simulation test case
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 7/13/2010
classdef Simulation_KnowledgeBase_Test < TestCase
    methods
        function this = Simulation_KnowledgeBase_Test(name)
            this = this@TestCase(name);
        end
        
        function testInitializationFromKnowledgeBase(~)
            dbConnectionParameters = config;
            database = edu.jiangnan.fmme.db.MySQLDatabase(dbConnectionParameters);
            kbWID = edu.jiangnan.fmme.cell.kb.KnowledgeBaseUtil.selectLatestKnowledgeBase(database);
            kb = edu.jiangnan.fmme.cell.kb.KnowledgeBase(database, kbWID);

            sim = edu.jiangnan.fmme.cell.sim.Simulation(kb.states, kb.processes);
            sim.initializeConstants(kb);
            sim.applyOptions('lengthSec', 10, 'verbosity', 0, 'seed', 1);
            fitter = edu.jiangnan.fmme.cell.sim.util.FitConstants(sim, struct('method', 'heuristic', 'verbosity', 0));
            fitter.run();
            sim.run();

            assertElementsAlmostEqual(...
                1 / (9.0 * 3600) * log(2), sim.state('MetabolicReaction').growth(:,:,1),...
                'relative', 35e-2);
            
            database.close();
        end
    end
end

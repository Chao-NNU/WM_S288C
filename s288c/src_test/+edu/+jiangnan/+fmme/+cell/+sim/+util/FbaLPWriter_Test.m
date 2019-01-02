%FbaLPWriter test
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updaetd 2/3/2011
classdef FbaLPWriter_Test < TestCase
    methods
        function this = FbaLPWriter_Test(name)
            this = this@TestCase(name);
        end
    end
    
    methods
        %checks for run time errors
        function testWriter(~)
            sim = edu.jiangnan.fmme.cell.sim.SimulationFixture.load([], true);
            edu.jiangnan.fmme.cell.sim.util.FbaLPWriter.write(sim, 'output/runSmallTests');
        end
    end
end
%ProteinComplex_Test
% Protein complex test class.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/4/2011
classdef ProteinComplex_Test < edu.jiangnan.fmme.cell.sim.CellStateTestCase
    methods
        function this = ProteinComplex_Test(name)
            this = this@edu.jiangnan.fmme.cell.sim.CellStateTestCase(name);
        end
    end
    
    methods
        function testMolecularWeights(this)
            import edu.jiangnan.fmme.util.ConstantUtil;
            
            s = this.state;
            sim = edu.jiangnan.fmme.cell.sim.SimulationFixture.load([], true);
            g = sim.gene;
            m = sim.state('Metabolite');
            r = sim.state('Rna');
            pm = sim.state('ProteinMonomer');
            
            rnaLens = sum(s.proteinComplexComposition(setdiff(1:end, g.mRNAIndexs), :, :), 3)' * ...
                max(0, r.lengths(r.matureIndexs(setdiff(1:end, r.matureMRNAIndexs))) - 1);
            monLens = sum(s.proteinComplexComposition(g.mRNAIndexs, :, :), 3)' * ...
                max(0, pm.lengths(pm.matureIndexs) - 1);
            
            assertElementsAlmostEqual(s.molecularWeights, ...
                s.baseCounts * m.molecularWeights ...
                - (ConstantUtil.elements.O +     ConstantUtil.elements.H) * repmat(rnaLens, 6, 1) ...
                - (ConstantUtil.elements.O + 2 * ConstantUtil.elements.H) * repmat(monLens, 6, 1), ...
                'relative', 1e-8, 0);
        end
    end
end
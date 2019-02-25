%ProteinComplex
%- nascent
%- mature
%- bound
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 1/5/2011
classdef ProteinComplex < edu.jiangnan.fmme.cell.sim.MoleculeCountState
    
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'molecularWeights'
            'baseCounts'
            'lengths'
            'halfLives'
            'compartments'
            'proteinComplexComposition'
            'minimumAverageExpression'
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'counts'
            };
        dependentStateNames     = {}; %names of properties which can be calculated from the simulation's state
    end
    
    %indices
    properties
        nascentIndexs           %index within complexs
        matureIndexs            %index within complexs
        inactivatedIndexs       %index within complexs
        boundIndexs             %index within complexs
        misfoldedIndexs         %index within complexs
        damagedIndexs           %index within complexs
        
        ribosome40SIndexs       %index within matureIndexs
        ribosome60SIndexs       %index within matureIndexs
        ribosome80SIndexs       %index within matureIndexs
        %translationFactorIndexs %index within matureIndexs
        ftsZGTPIndexs           %index within matureIndexs
        ftsZGDPIndexs           %index within matureIndexs
		rnaPolymeraseIndexs     %index within matureIndexs
        %rnaPolymerase2Indexs    %index within matureIndexs
		%rnaPolymerase3Indexs    %index within matureIndexs
        replisomeIndexs         %index within matureIndexs
		dnaPolymeraseIndexs    %index within matureIndexs
		%dnaPolymeraseAIndexs    %index within matureIndexs
        %dnaPolymeraseDIndexs    %index within matureIndexs
		%dnaPolymeraseEIndexs    %index within matureIndexs
		%dnaPolymeraseGIndexs    %index within matureIndexs
		%dnaPolymeraseZIndexs    %index within matureIndexs
		ORCIndexs				%index within matureIndexs
        %dnaAPolymerIndexs       %index within matureIndexs
    end
    
    properties
        proteinComplexComposition %protein complex composition (monomers X complexes X compartments)
        formationProcesses        %indices of proceses where each complex is formed
        minimumAverageExpression  %minimum average expression
    end
    
    %references to objects
    properties
        chromosome
        rnaPolymerase
        ribosome
        ftsZRing
    end
    
    %constructor
    methods
        function this = ProteinComplex(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.MoleculeCountState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.jiangnan.fmme.cell.sim.MoleculeCountState(simulation);            
            this.chromosome = simulation.state('Chromosome');
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.ribosome = simulation.state('Ribosome');
            this.ftsZRing = simulation.state('FtsZRing');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.MoleculeCountState(knowledgeBase, simulation);
            
            numComplexs = knowledgeBase.numProteinComplexs;
            
            this.nascentIndexs     = (1:numComplexs)';
            this.matureIndexs      = (1:numComplexs)' + this.nascentIndexs(end);
            this.inactivatedIndexs = (1:numComplexs)' + this.matureIndexs(end);
            this.boundIndexs       = (1:numComplexs)' + this.inactivatedIndexs(end);
            this.misfoldedIndexs   = (1:numComplexs)' + this.boundIndexs(end);
            this.damagedIndexs     = (1:numComplexs)' + this.misfoldedIndexs(end);
            
            this.wholeCellModelIDs = repmat({knowledgeBase.proteinComplexs.wholeCellModelID}', 6, 1);
            this.names = repmat({knowledgeBase.proteinComplexs.name}', 6, 1);
            this.baseCounts = [
                reshape([knowledgeBase.proteinComplexs.baseCount], [], numComplexs)';
                repmat(reshape([knowledgeBase.proteinComplexs.matureBaseCount], [], numComplexs)', 5, 1)];
            this.molecularWeights = [...
                knowledgeBase.proteinComplexs.molecularWeight  ...
                repmat([knowledgeBase.proteinComplexs.matureMolecularWeight], 1, 5)]';
            this.compartments = double(repmat(knowledgeBase.proteinComplexCompartments, 6, 1));
            this.halfLives = [...
                knowledgeBase.proteinComplexs.halfLife  ...
                repmat([knowledgeBase.proteinComplexs.matureHalfLife], 1, 3) ...
                zeros(1, knowledgeBase.numProteinComplexs) ...
                zeros(1, knowledgeBase.numProteinComplexs)]';
            
            this.proteinComplexComposition = knowledgeBase.proteinComplexAllRNAComposition;
            this.proteinComplexComposition([knowledgeBase.mRNAGenes.idx], :, :) = knowledgeBase.proteinComplexAllMonomerComposition;
            
            tmp = [knowledgeBase.proteinComplexs.complexFormationProcess];
            this.formationProcesses = repmat([tmp.idx]', 6, 1);
            
            %this.halfLives(~ismember(this.compartments, [this.compartment.cytosolIndexs; this.compartment.terminalOrganelleCytosolIndexs])) = Inf;
            this.halfLives(~ismember(this.compartments, [this.compartment.cytosolIndexs; this.compartment.mitochondriaIndexs])) = Inf;
			
            this.ribosome40SIndexs = this.getIndexs('RIBOSOME_40S');
            this.ribosome60SIndexs = this.getIndexs('RIBOSOME_60S');
            this.ribosome80SIndexs = this.getIndexs('RIBOSOME_80S');
            %{
			this.translationFactorIndexs = this.getIndexs({
				'YJR007W_YPL237W_YER025W_TRIMER'
				'YBR079C_YDR429C_YLR192C_YMR146C_YMR309C_YOR361C_6MER'
                });
			%}	
            this.ftsZGTPIndexs = this.getIndexs('YDL126C_9MER_GTP');
            this.ftsZGDPIndexs = this.getIndexs('YDL126C_9MER_GDP');
            this.rnaPolymeraseIndexs = this.getIndexs({
                'RNA_POLYMERASE_II'
				'RNA_POLYMERASE_III'
                });
            this.replisomeIndexs = this.getIndexs({
                'YLR274W_YBR202W_YPR019W_YBL023C_YEL032W_YGL201C_6MER'
                'YOR217W_YJR068W_YNL290W_YOL094C_YBR087W_PENTAMER'
                %'DNA_POLYMERASE_CORE_BETA_CLAMP_GAMMA_COMPLEX'
                %'DNA_POLYMERASE_CORE_BETA_CLAMP_PRIMASE'
                });
            this.dnaPolymeraseIndexs = this.getIndexs({
				'DNA_POLYMERASE_A'
				'DNA_POLYMERASE_D'
				'DNA_POLYMERASE_E'
				'DNA_POLYMERASE_Z'                
                });
			this.ORCIndexs = this.getIndexs({
				'YML065W_1MER_ATP'      										%Orc1-ATP 1mer
                'YML065W_YBR060C_2MER_ATP_ADP' 									%Orc 2mer-(1)ATP-(1)ADP
                'YML065W_YBR060C_2MER_2ATP'      								%Orc-ATP 2mer
                'YML065W_YBR060C_YLL004W_3MER_2ATP_ADP' 						%Orc 3mer-(2)ATP-(1)ADP
                'YML065W_YBR060C_YLL004W_3MER_3ATP'      						%Orc-ATP 3mer
                'YML065W_YBR060C_YLL004W_YPR162C_4MER_3ATP_ADP' 				%Orc 4mer-(3)ATP-(1)ADP
                'YML065W_YBR060C_YLL004W_YPR162C_4MER_4ATP'      				%Orc-ATP 4mer
                'YML065W_YBR060C_YLL004W_YPR162C_YNL261W_5MER_4ATP_ADP' 		%Orc 5mer-(4)ATP-(1)ADP
                'YML065W_YBR060C_YLL004W_YPR162C_YNL261W_5MER_5ATP'      		%Orc-ATP 5mer
                'YML065W_YBR060C_YLL004W_YPR162C_YNL261W_YHR118C_6MER_5ATP_ADP' %Orc 6mer-(5)ATP-(1)ADP
                'YML065W_YBR060C_YLL004W_YPR162C_YNL261W_YHR118C_6MER_6ATP'     %Orc-ATP 6mer
				});
        end
    end
            
    methods
        function notUpdatingProteins = updateExternalState(this, deltaProteins, proteinIsDegraded)
            c = this.chromosome;
			
				notUpdatingProteins = zeros(size(deltaProteins));
				deltaFreeProteins = deltaProteins(this.matureIndexs, this.compartment.cytosolIndexs);
				deltaBoundProteins = deltaProteins(this.boundIndexs, this.compartment.cytosolIndexs);

				%update ribosome state
				this.ribosome.releaseRibosome(-deltaBoundProteins(this.ribosome80SIndexs), 0);
				deltaBoundProteins(this.ribosome80SIndexs) = 0;
            
				%update bound translation elongation factors
				%deltaBoundProteins(this.translationFactorIndexs) = 0;
            
				%update bound FtsZ state
				notUpdatingProteins(this.boundIndexs([this.ftsZGTPIndexs; this.ftsZGDPIndexs]), this.compartment.cytosolIndexs) = ...
					this.ftsZRing.releaseFtsZ(-deltaBoundProteins([this.ftsZGTPIndexs; this.ftsZGDPIndexs]));
				deltaBoundProteins(this.ftsZGTPIndexs) = 0;
				deltaBoundProteins(this.ftsZGDPIndexs) = 0;
				
				%prevent changes to bound DnaA complexes to that replication
				%initiation not disrupted
				if any(deltaBoundProteins(this.ORCIndexs) < 0)
					warning('WholeCell:warning', 'Origin recognition complex not decayed');
					notUpdatingProteins(this.boundIndexs(this.ORCIndexs), this.compartment.cytosolIndexs) = ...
						-deltaBoundProteins(this.ORCIndexs);
					deltaBoundProteins(this.ORCIndexs) = 0;
				end
				
				%prevent changes to bound DNA polymerase so that replication
				%not disrupted; issue warning
				if any(deltaBoundProteins(this.replisomeIndexs) < 0)
					warning('WholeCell:warning', 'DNA polymerase not decayed');
					notUpdatingProteins(this.boundIndexs(this.replisomeIndexs), this.compartment.cytosolIndexs) = ...
						-deltaBoundProteins(this.replisomeIndexs);
					deltaBoundProteins(this.replisomeIndexs) = 0;
				end
				
				%update chromosomally bound proteins
				idxs = find(deltaBoundProteins < 0);
				
			for i = 16				
				[posStrnds{i}, proteins{i}] = find(c.complexBoundSites{i});
				chrReleasePosStrnds{i} = zeros(0, 2);
				for j = 1:numel(idxs)
					if sum(proteins{i} == idxs(j)) < -deltaBoundProteins(idxs(j))
						%throw(MException('ProteinComplex:error', 'Error updating external state'))
						warning('Error updating external state');
					end
					
					chrReleasePosStrnds{i} = [
						chrReleasePosStrnds{i};
						this.randStream.randomlySelectNRows(posStrnds{i}(proteins{i} == idxs(j), :), -deltaBoundProteins(idxs(j)))
						]; %#ok<AGROW>
				end
				
				[~, releasedComplexs{i}] = c.setRegionProteinUnbound(i,...
					chrReleasePosStrnds{i}, 1, [], idxs, ...
					false, false, false, proteinIsDegraded);
									
				if ~isequal(-deltaBoundProteins(idxs), releasedComplexs{i})
					%throw(MException('ProteinComplex:error', 'Chromosomally bound proteins impropely released'));
					warning('Chromosomally bound proteins impropely released');
				end
			end				
				
				%update free RNA polymerase
			if deltaFreeProteins(this.rnaPolymeraseIndexs(2)) && proteinIsDegraded
				this.rnaPolymerase.degradeFreePolymerase(-deltaFreeProteins(this.rnaPolymeraseIndexs(2)));
			end
		end
    end
    
    %helper methods
    methods
        function value = getIndexs(this, wholeCellModelIDs)
            [~, value] = ismember(wholeCellModelIDs, this.wholeCellModelIDs(this.matureIndexs));
        end
    end
end

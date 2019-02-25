%tRNA Aminoacylation
%
% @wholeCellModelID Process_tRNAAminoacylation
% @name             tRNA Aminoacylation
% @description
%   Biology
%   ===============
%   In biological systems tRNAs serve as mediators between the ribosome and
%   the amino acids which forms peptide polymers. In this process we
%   simulate the conjugation of amino acids to the tRNAs which deliver them
%   to the ribosome. Additionally we simulate the aminoacylation of the
%   tmRNA which similarly delivers the amino acid alanine to stalled
%   ribosomes.
%
%   Knowledge Base
%   ===============
%   As of 8/11/2010 the M. genitalium knowledge base contains 39 tRNA
%   aminoacylation reactions involving
%   - 37 aminoacylation reactions, each assuming a cost of 1 ATP per
%     aminoacylation
%   - 2 transfer reactions
%   - 36 tRNAs and 1 tmRNA
%   - 21 enzymes
%   - 20 amino acid substrates and 10 additional substrates
%
%   The reaction kinetics stored in the knowledge base were compiled
%   several sources including SABIO-RK [PUB_0100].
%
%   Note: Unlike knowledge base representation, for transfer reactions
%   reactionStoichiometryMatrix doesn't include the already conjugated
%   amino acid on the left-hand-side or the the ultimate conjugated amino
%   acid on the right-hand-side.
%
%   Representation
%   ===============
%   substrates represents the counts of free metabolites available for tRNA
%   aminoacylation. enzymes represents the counts of proteins available to
%   catalyze tRNA aminoacylation reactions. freeRNAs and aminoacylatedRNAs
%   represent the counts of free, unaminoacylated and aminocylated RNAs. We do
%   not represent intermediate states in the aminoacylation of tRNAs.  The
%   molecular weights of the unaminoacylated and aminocylated tRNAs are computed
%   by the knowledge base RNA classes.
%
%   reactionStoichiometryMatrix, reactionModificationMatrix,
%   reactionCatalysisMatrix, and enzymeBounds represent the tRNA aminoacylation
%   reactions. reactionStoichiometryMatrix represents the free metabolites
%   required for  each reaction. reactionModificationMatrix represents the tRNA
%   aminoacylated by each reaction. reactionCatalysis represents the enzyme
%   required to aminoacylate each tRNA. enzymeBounds represents the kcat of the
%   catalyzing enzyme of each reaction. Note: reactionStoichiometryMatrix here
%   is slightly modified from that computed by the super class: amino acids
%   conjugated to tRNAs before and after transfer reactions have been removed
%   from the left- and righ-hand-side of reactionStoichiometryMatrix.
%
%   Initialization
%   ===============
%   tRNAs are all initialized to the aminoacylated state.
%
%   Simulation
%   ===============
%   Uses greedy algorithm to simulate tRNA aminacylation (complexation of
%   tRNAs with the specific amino acids that aminoacylate them).
%   1. Deterministically activate tRNAs up to the minimum of free tRNAs and
%      amino acids, proportional to free tRNAs
%   2. Stoichastically activate residual tRNAs using residual amino acids
%      with probabilities proportional to remaining free tRNAs
%
%      While(true)
%        1. Calculate numbers of protein monomers that can be modified based on
%           substrate, enzyme, and unmodified protein monomer availability and
%           kinetics.
%        2. Randomly select protein monomer to modify weighted by limits
%           calculated in step (1).
%        3. Update substrates, enzymes, unmodified protein monomers
%        4. Repeat until insufficient resources to further modify protein
%           monomers
%      End
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanfod.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 8/11/2010

classdef tRNAAminoacylation < edu.jiangnan.fmme.cell.sim.ReactionProcess
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
		fixedConstantNames__       = {   %names of fixed constant properties
			'monomerTRNACounts';
			};
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'freeRNAs';
            'aminoacylatedRNAs'};
    end

    %IDs, names, and local indices
    properties
        substrateIndexs_aminoAcids                   %indices within substrates of amino acids
        substrateIndexs_glutamate                    %index within substrates of glutamate
        substrateIndexs_glutamine                    %index within substrates of glutamine
        substrateIndexs_methionine                   %index within substrates of methionine
        substrateIndexs_fmethionine                  %index within substrates of formylmethionine
        substrateIndexs_atp                          %index within substrates of ATP
        substrateIndexs_amp                          %index within substrates of AMP
        substrateIndexs_adp                          %index within substrates of ADP
        substrateIndexs_diphosphate                  %index within substrates of diphosphate
        substrateIndexs_phosphate                    %index within substrates of inorganic phosphate
        substrateIndexs_water                        %index within substrates of water
        substrateIndexs_hydrogen                     %index within substrates of hydrogen
        substrateIndexs_fthf10                       %index within substrates of 10-forymltetrahydrofolate
        substrateIndexs_thf                          %index within substrates of tetrahydrofolate

        enzymeIndexs_tRNASynthetases                 %indices within enzymes of tRNA synthetases
        enzymeIndexs_tRNATransferases                %indices within enzymes of tRNA transferases
        enzymeIndexs_tRNAGlutamylSynthetase          %index within enzymes of tRNA glutamyl synthetase
        enzymeIndexs_tRNAMethionylSynthetase         %index within enzymes of tRNA methionyl synthetase
        enzymeIndexs_tRNAGlutamylAmidotransferase    %index within enzymes of tRNA glutamyl amidotransferase
        enzymeIndexs_tRNAMethionylFormyltransferase  %index within enzymes of tRNA methionyl formyltransferase
        enzymeIndexs_tRNAs                           %index within enzymes of tRNAs
        %enzymeIndexs_tmRNA                           %index within enzymes of tmRNA

        reactionIndexs_aminoacylation                %indices within reactions of tRNA aminoacylation reactions
        reactionIndexs_transfer                      %indices within reactions of tRNA transfer reactions
        reactionIndexs_glutamylamidotransfer         %index within reactions of tRNA glutamylamidotransfer reaction
        reactionIndexs_methionylformyltransfer       %index within reactions of tRNA methionylformyltransfer reaction

        compartmentWholeCellModelIDs                 %whole cell model ids of compartments
        compartmentIndexs_cytosol                    %index within compartments of cytosol
		compartmentIndexs_mitochondrion              %index within compartments of mitochondrial

        freeRNAWholeCellModelIDs                     %whole cell model ids of free RNAs
        aminoacylatedRNAWholeCellModelIDs            %whole cell model ids of aminoacylated RNAs

        aminoacylatedRNAGlobalIndexs                 %indices of aminoacylated RNAs within simulation.matureRNAIndexs
        aminoacylatedRNAtRNAIndexs                   %indices of aminoacylated tRNAs within aminoacylatedRNAGlobalIndexs
        
        speciesIndexs_enzymes                        %indexs of enzyme species within speciesReactantByproductMatrix, speciesReactantMatrix
    end

    %fixed biological constants
    properties
        monomerTRNACounts              %numbers of each tRNA species required to translate each protein monomer species    
        
        speciesReactantByproductMatrix %stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; free RNAs]
        speciesReactantMatrix          %reactant stoichiometry of susbtrates, enzymes, RNA in modifications necessary to mature each modifying RNA: [RNA modifications] X [subtrates; enzymes; free RNAs]
    end

    %global state (stored locally for convenience)
    properties
        freeRNAs          %counts of free RNAs
        aminoacylatedRNAs %counts of aminoacylated RNAs
    end
    
    %dependent global state (implemented as dependent property for
    %convenience)
    properties (Dependent)
        tRNASynthetases   %counts of tRNA synthetases
        tRNATransferases  %counts of tRNA transferases
    end

    %constructor
    methods
        function this = tRNAAminoacylation(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.ReactionProcess(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            %import classes
            import edu.jiangnan.fmme.util.ComputationUtil;

            %call super class method to get knowledge base reactions
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.ReactionProcess(knowledgeBase, simulation, varargin{:});

            %compartments
            this.compartmentWholeCellModelIDs = this.compartment.wholeCellModelIDs(this.compartment.cytosolIndexs); 
            this.compartmentIndexs_cytosol = 1;
			this.compartmentIndexs_mitochondrion = 4;%add by manual

            %create representation for GLN which super class leaves off
            %because it cancels out in the knowledge base
            this.substrateWholeCellModelIDs  = [this.substrateWholeCellModelIDs; 'GLN_L'];
            this.reactionStoichiometryMatrix = [this.reactionStoichiometryMatrix; zeros(1, size(this.reactionStoichiometryMatrix, 2))];
            this.reactionCoenzymeMatrix      = [this.reactionCoenzymeMatrix zeros(size(this.reactionCoenzymeMatrix, 1), 1)];

            %substrate indices
            this.substrateIndexs_aminoAcids  = this.substrateIndexs({'ALA_L'; 'ARG_L'; 'ASN_L'; 'ASP_L'; 'CYS_L'; 'GLN_L'; 'GLU_L'; 'GLY_L'; 'HIS_L'; 'ILE_L'; 'LEU_L'; 'LYS_L'; 'MET_L'; 'PHE_L'; 'PRO_L'; 'SER_L'; 'THR_L'; 'TRP_L'; 'TYR_L'; 'VAL_L'});
            this.substrateIndexs_glutamate   = this.substrateIndexs({'GLU_L'});
            this.substrateIndexs_glutamine   = this.substrateIndexs({'GLN_L'});
            this.substrateIndexs_methionine  = this.substrateIndexs({'MET_L'});
            this.substrateIndexs_fmethionine = this.substrateIndexs({'FMET'});
            this.substrateIndexs_atp         = this.substrateIndexs({'ATP'});
            this.substrateIndexs_amp         = this.substrateIndexs({'AMP'});
            this.substrateIndexs_adp         = this.substrateIndexs({'ADP'});
            this.substrateIndexs_diphosphate = this.substrateIndexs({'PPI'});
            this.substrateIndexs_phosphate   = this.substrateIndexs({'PI'});
            this.substrateIndexs_water       = this.substrateIndexs({'H2O'});
            this.substrateIndexs_hydrogen    = this.substrateIndexs({'H'});
            this.substrateIndexs_fthf10      = this.substrateIndexs({'10FTHF'});
            this.substrateIndexs_thf         = this.substrateIndexs({'THF'});

            %include tRNAs as enzymes
            this.enzymeWholeCellModelIDs = [
                this.enzymeWholeCellModelIDs;
                %this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureTRNAIndexs));
                %'MG_0004'];
				this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureTRNAIndexs))];
				
            %enzyme indices
			this.enzymeIndexs_tRNASynthetases = this.enzymeIndexs({
				'YOR335C_MONOMER';		%alanyl-tRNA synthetase
				'YBL076C_MONOMER';		%Isoleucyl-tRNA synthetase, cytoplasmic
				'YPL040C_MONOMER';		%Isoleucyl-tRNA synthetase, mitochondrial
				'YDR023W_MONOMER';		%Seryl-tRNA synthetase, cytoplasmic
				'YHR011W_MONOMER';		%Seryl-tRNA synthetase, mitochondrial
				'YIL078W_MONOMER';		%Threonyl-tRNA synthetase, cytoplasmic
				'YKL194C_MONOMER';		%Threonyl-tRNA synthetase, mitochondrial
				'YNL247W_MONOMER';		%Cysteinyl-tRNA synthetase
				'YHR020W_MONOMER';		%Prolyl-tRNA synthetase, cytoplasmic
				'YER087W_MONOMER';		%Prolyl-tRNA synthetase, mitochondrial
				'YGR264C_MONOMER';		%Methionyl-tRNA synthetase, cytoplasmic
				'YGR171C_MONOMER';		%Methionyl-tRNA synthetase, mitochondrial
				'YLL018C_MONOMER';		%Aspartyl-tRNA synthetase, cytoplasmic
				'YPL104W_MONOMER';		%Aspartyl-tRNA synthetase, mitochondrial
				'YLR060W_YFL022C_DIMER';%Phenylalanyl-tRNA synthetase, cytoplasmic
				'YPR047W_MONOMER';		%Phenylalanyl-tRNA synthetase, mitochondrial
				'YDR341C_MONOMER';		%Arginyl-tRNA synthetase, cytoplasmic
				'YHR091C_MONOMER';		%Arginyl-tRNA synthetase, mitochondrial
				'YPR081C_MONOMER';		%Glycyl-tRNA synthetase, cytoplasmic
				'YBR121C_MONOMER';		%Glycyl-tRNA synthetase, mitochondrial
				'YPL160W_MONOMER';		%Leucyl-tRNA synthetase, cytoplasmic
				'YLR382C_MONOMER';		%Leucyl-tRNA synthetase, mitochondrial
				'YDR037W_MONOMER';		%Lysyl-tRNA synthetase, cytoplasmic
				'YNL073W_MONOMER';		%Lysyl-tRNA synthetase, mitochondrial
				'YGL245W_MONOMER';		%Glutamyl-tRNA synthetase, cytoplasmic
				'YOL033W_MONOMER';		%Glutamyl-tRNA synthetase, mitochondrial
				'YGR185C_MONOMER';		%Tyrosyl-tRNA synthetase, cytoplasmic
				'YPL097W_MONOMER';		%Tyrosyl-tRNA synthetase, mitochondrial
				'YOL097C_MONOMER';		%Tryptophanyl-tRNA synthetase, cytoplasmic
				'YDR268W_MONOMER';		%Tryptophanyl-tRNA synthetase, mitochondrial
				'YGR094W_MONOMER';		%Valyl-tRNA synthetase, mitochondrial
				'YHR019C_MONOMER';		%Asparaginyl-tRNA synthetase, cytoplasmic
				'YCR024C_MONOMER';		%Asparaginyl-tRNA synthetase, mitochondrial
				'YPR033C_MONOMER';		%Histidyl-tRNA synthetase, mitochondrial
				'YOR168W_MONOMER'});	%Glutaminyl-tRNA synthetase

				this.enzymeIndexs_tRNATransferases = this.enzymeIndexs({
				'YMR293C_YBL080C_YGR102C_6MER';	    %glutamyl-tRNA(Gln) amidotransferase complex
				'YBL013W_MONOMER'});				%Formyl-Methionyl-tRNA Transformylase
			
			this.enzymeIndexs_tRNAGlutamylSynthetase         = this.enzymeIndexs({'YGL245W_MONOMER'; 'YOL033W_MONOMER'});  %glutamyl-tRNA synthetase
			this.enzymeIndexs_tRNAMethionylSynthetase        = this.enzymeIndexs({'YGR264C_MONOMER'; 'YGR171C_MONOMER'});  %methionyl-tRNA synthetase
			this.enzymeIndexs_tRNAGlutamylAmidotransferase   = this.enzymeIndexs({'YMR293C_YBL080C_YGR102C_6MER'}); %glutamyl-tRNA(Gln) and/or aspartyl-tRNA(Asn)
			
			%amidotransferase
			this.enzymeIndexs_tRNAMethionylFormyltransferase = this.enzymeIndexs({'YBL013W_MONOMER'});        %methionyl-tRNA formyltransferase
            this.enzymeIndexs_tRNAs                          = this.enzymeIndexs(this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.rna.matureTRNAIndexs)));
			
            %RNAs
            this.reactionModificationMatrix(isnan(this.reactionModificationMatrix)) = 1;
            
			aaGenes = any(this.reactionModificationMatrix, 1)';
			aaRNAs  = ComputationUtil.invertCompositionMatrix(this.rna.matureRNAGeneComposition) * aaGenes;
			
            this.aminoacylatedRNAGlobalIndexs  = find(aaRNAs);
            this.aminoacylatedRNAtRNAIndexs    = intersect(this.aminoacylatedRNAGlobalIndexs, this.rna.matureTRNAIndexs);

            this.freeRNAWholeCellModelIDs          = this.rna.wholeCellModelIDs(this.rna.matureIndexs(this.aminoacylatedRNAGlobalIndexs));
            this.aminoacylatedRNAWholeCellModelIDs = this.rna.wholeCellModelIDs(this.rna.aminoacylatedIndexs(this.aminoacylatedRNAGlobalIndexs));

            this.reactionModificationMatrix = this.reactionModificationMatrix(:, aaGenes);

            %reaction indices
            this.reactionIndexs_aminoacylation          = find(strcmp(this.reactionTypes, 'aminoacylation'));
            this.reactionIndexs_transfer                = find(strcmp(this.reactionTypes, 'transfer'));
            %this.reactionIndexs_glutamylamidotransfer   = this.reactionIndexs({'MG502_Amidotransferase'});
			this.reactionIndexs_glutamylamidotransfer   = this.reactionIndexs({
			'tQ(CUG)M_Amidotransferase';
			'tQ(UUG)B_Amidotransferase';
			'tQ(UUG)C_Amidotransferase'; 
			'tQ(UUG)D1_Amidotransferase';
			'tQ(UUG)D2_Amidotransferase';
			'tQ(UUG)D3_Amidotransferase';
			'tQ(UUG)E1_Amidotransferase'; 
			'tQ(UUG)E2_Amidotransferase';
			'tQ(UUG)H_Amidotransferase'; 
			'tQ(UUG)L_Amidotransferase';
			'tQ(UUG)Q_Amidotransferase'});
            %this.reactionIndexs_methionylformyltransfer = this.reactionIndexs({'MG488_Formyltransferase'});
			this.reactionIndexs_methionylformyltransfer = this.reactionIndexs({'tM(CAU)J1_Formyltransferase'});

            %for transfer reactions unlike knowledge base representation,
            %formulate reactionStoichiometryMatrix so that left-hand-side
            %doesn't include the already conjugated amino acid, and the
            %right-hand-side doesn't include the ultimate conjugated amino
            %acid
            this.reactionStoichiometryMatrix(this.substrateIndexs_glutamine,   this.reactionIndexs_glutamylamidotransfer)   = -1;
            this.reactionStoichiometryMatrix(this.substrateIndexs_glutamate,   this.reactionIndexs_glutamylamidotransfer)   =  1;
            this.reactionStoichiometryMatrix(this.substrateIndexs_methionine,  this.reactionIndexs_methionylformyltransfer) =  0;
            this.reactionStoichiometryMatrix(this.substrateIndexs_fmethionine, this.reactionIndexs_methionylformyltransfer) =  0;

            %re-run super class method to update substrate and enzyme
            %mappings to include GLN, the tRNAs, and the tmRNA
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.Process(knowledgeBase, simulation, varargin{:});
            
            this.initializeSpeciesNetwork();

            %protein monomers -- needed for fitting
            this.monomerTRNACounts = knowledgeBase.proteinMonomerTRNACounts;
        end
        
        %initialize reaction stoichiomety to be used in evolveState:
        %  [tRNA aminoacylations and transfers] X [subtrates; enzymes; free tRNAs]
        function initializeSpeciesNetwork(this)
            % formulate network: [tRNA aminoacylations and transfers] X [subtrates; enzymes; free tRNAs]
            numMetabolites = size(this.reactionStoichiometryMatrix, 1);
            numEnzymes     = size(this.reactionCatalysisMatrix, 2);
            numRNAs        = size(this.reactionModificationMatrix, 2);
            
            this.speciesReactantByproductMatrix = [this.reactionModificationMatrix' * [ ...
                -this.reactionStoichiometryMatrix' ...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:, 2) * this.stepSizeSec, 1, numEnzymes)] ...
                eye(numRNAs)];
            this.speciesReactantMatrix = [this.reactionModificationMatrix' * max(0, [ ...
                -this.reactionStoichiometryMatrix' ...
                this.reactionCatalysisMatrix ./ repmat(this.enzymeBounds(:, 2) * this.stepSizeSec, 1, numEnzymes)]) ...
                eye(numRNAs)];
            
            this.speciesIndexs_enzymes = (1:numEnzymes)' + numMetabolites;
        end

        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.jiangnan.fmme.cell.sim.ReactionProcess();

            this.freeRNAs          = this.rna.counts(this.rna.matureIndexs(this.aminoacylatedRNAGlobalIndexs), this.compartment.cytosolIndexs, :);
            this.aminoacylatedRNAs = this.rna.counts(this.rna.aminoacylatedIndexs(this.aminoacylatedRNAGlobalIndexs), this.compartment.cytosolIndexs, :);
        end

        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.jiangnan.fmme.cell.sim.ReactionProcess();
            this.rna.counts(this.rna.matureIndexs(       this.aminoacylatedRNAGlobalIndexs), this.compartment.cytosolIndexs, :) = this.freeRNAs;
            this.rna.counts(this.rna.aminoacylatedIndexs(this.aminoacylatedRNAGlobalIndexs), this.compartment.cytosolIndexs, :) = this.aminoacylatedRNAs;
        end
    end

    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.jiangnan.fmme.cell.sim.ReactionProcess(numTimePoints);

            this.freeRNAs          = zeros(length(this.aminoacylatedRNAGlobalIndexs), 1, numTimePoints);
            this.aminoacylatedRNAs = zeros(length(this.aminoacylatedRNAGlobalIndexs), 1, numTimePoints);
        end
    end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% substrate and byproducts
            %aminoacylations of each t(m)RNA
            [~, idxs] = ismember(this.aminoacylatedRNAtRNAIndexs, this.aminoacylatedRNAGlobalIndexs);
            aminoacylations = zeros(size(this.aminoacylatedRNAGlobalIndexs));
			
            aminoacylations(idxs) = aminoacylations(idxs) + ...
                this.monomerTRNACounts' * states.monomerProductions;
            
            %t(m)RNA aminoacylation reactions
            reactions = this.reactionModificationMatrix * aminoacylations;
            
            %metabolic load of t(m)RNA aminoacylation reactions
            bmProd = max(0, -this.reactionStoichiometryMatrix) * reactions;
            byProd = max(0,  this.reactionStoichiometryMatrix) * reactions;
            
            %% enzymes
            %initialize
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %tRNAs
            minEnzExp(this.enzymeIndexs_tRNAs) = 1.95 * this.monomerTRNACounts' * states.monomerProductions0;
            
            %tmRNA
            %minEnzExp(this.enzymeIndexs_tmRNA) = 0;

            %tRNA synthetases, transferases
            [~, idxs] = ismember(this.aminoacylatedRNAtRNAIndexs, this.aminoacylatedRNAGlobalIndexs);
            aminoacylations = states.rnaProductions0(this.aminoacylatedRNAGlobalIndexs);
            aminoacylations(idxs)                = max(aminoacylations(idxs),                minEnzExp(this.enzymeIndexs_tRNAs));
            %aminoacylations(setdiff(1:end,idxs)) = max(aminoacylations(setdiff(1:end,idxs)), minEnzExp(this.enzymeIndexs_tmRNA));            
            
            minEnzExp(setdiff(1:end, [this.enzymeIndexs_tRNAs])) = ...
				+ minEnzExp(setdiff(1:end, [this.enzymeIndexs_tRNAs])) ...
				+ 3 * this.reactionCatalysisMatrix' * ((this.reactionModificationMatrix * aminoacylations) ./ this.enzymeBounds(:, 2));
			
			%{
			minEnzExp(setdiff(1:end, [this.enzymeIndexs_tRNAs; this.enzymeIndexs_tmRNA])) = ...
                + minEnzExp(setdiff(1:end, [this.enzymeIndexs_tRNAs; this.enzymeIndexs_tmRNA])) ...
                + 3 * this.reactionCatalysisMatrix' * ((this.reactionModificationMatrix * aminoacylations) ./ this.enzymeBounds(:, 2));
			%}
        end

        %initialization: RNAs initialized to mature state by simulation
        %initializeState method. This class initializes RNAs to
        %aminoacylated state.
        function initializeState(this)
            totalRNAs = this.aminoacylatedRNAs + this.freeRNAs;
            this.aminoacylatedRNAs = ceil(2 / 3 * totalRNAs);
            this.freeRNAs = totalRNAs - this.aminoacylatedRNAs;
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = max(0, -this.reactionStoichiometryMatrix) * min(...
                this.reactionModificationMatrix * (this.freeRNAs + this.aminoacylatedRNAs + 1), ...
                this.reactionCatalysisMatrix * this.enzymes(1:numel(this.speciesIndexs_enzymes)) * this.stepSizeSec);
        end
        
        %simulation
        function evolveState(this)
            %terminate early if no free RNAs
            if ~any(this.freeRNAs)
                return;
            end
            
            numRNAs = numel(this.freeRNAs);
            
            species = [
                this.substrates;
                this.enzymes(1:numel(this.speciesIndexs_enzymes));
                this.freeRNAs]';
            
            reactionFluxes = zeros(size(this.speciesReactantByproductMatrix, 1), 1);
            anyFlux = false;
            
            limits = species(ones(numRNAs, 1), :) ./ this.speciesReactantMatrix;
            limits(:, this.substrateIndexs_water) = NaN;
            limits(:, this.substrateIndexs_hydrogen) = NaN;
            reactionLimits = min(limits, [], 2)';
            reactionLimits(isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 0) = 0;
            isReactionInactive = reactionLimits <= 0;
            
            while true
                %compute maximum number of each species that can be modified
                limits = species(ones(numRNAs, 1), :) ./ max(0, this.speciesReactantByproductMatrix);                
                limits(:, this.substrateIndexs_water) = NaN;
                limits(:, this.substrateIndexs_hydrogen) = NaN;
                reactionLimits = min(...
                    this.randStream.stochasticRound(min(limits(:, this.speciesIndexs_enzymes), [], 2)), ...
                    min(limits(:, [1:this.speciesIndexs_enzymes(1)-1  this.speciesIndexs_enzymes(end)+1:end]), [], 2))';
                reactionLimits(isReactionInactive | isinf(reactionLimits) | isnan(reactionLimits) | reactionLimits < 1) = 0;
                
                %stop if no more tRNA can be aminoacylated
                if ~any(reactionLimits)
                    break;
                end
                anyFlux = true;
                               
                edges = min([0 cumsum(reactionLimits(:)' / sum(reactionLimits))], 1);
                nRxns = min(reactionLimits(reactionLimits > 0));
                if nRxns <= 1
                    %pick substrate, reaction, enzymes
                    [~, selectedReaction] = histc(rand(this.randStream, 1, 1), edges);
                    reactionFluxes(selectedReaction) = ...
                        + reactionFluxes(selectedReaction) ...
                        + 1;
                    
                    %decrement metabolites, enzymes, unmodified substrates; increment modified substrates
                    species = species - this.speciesReactantByproductMatrix(selectedReaction, :);
                else
                    %pick substrate, reaction, enzymes
                    edges(end) = 1.1;
                    selectedReactions = histc(rand(this.randStream, nRxns, 1), edges);
                    selectedReactions = selectedReactions(1:end-1);
                    selectedReactions = selectedReactions(:);
                    reactionFluxes = ...
                        + reactionFluxes ...
                        + selectedReactions;
                    
                    %decrement metabolites, enzymes, unmodified substrates; increment modified substrates
                    species = species - selectedReactions' * this.speciesReactantByproductMatrix;
                end
            end
            
            %stop if no reaction can proceed
            if ~anyFlux
                return;
            end
            
            % update substrates
            this.substrates = this.substrates + ...
                this.reactionStoichiometryMatrix * this.reactionModificationMatrix * reactionFluxes;
            
            % update RNAs
            this.freeRNAs          = this.freeRNAs          - reactionFluxes;
            this.aminoacylatedRNAs = this.aminoacylatedRNAs + reactionFluxes;
        end
    end

    %get methods of dependent global state
    methods
        function value = get.tRNASynthetases(this)
            value = this.enzymes(this.enzymeIndexs_tRNASynthetases,:,:);
        end

        function value = get.tRNATransferases(this)
            value = this.enzymes(this.enzymeIndexs_tRNATransferases,:,:);
        end
    end

    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.freeRNAs, 3) == 1
                value = this.getDryWeight@edu.jiangnan.fmme.cell.sim.ReactionProcess() + (...
                    this.rna.molecularWeights(this.rna.matureIndexs(this.aminoacylatedRNAGlobalIndexs))'        * this.freeRNAs         + ...
                    this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.aminoacylatedRNAGlobalIndexs))' * this.aminoacylatedRNAs) ...
                    / edu.jiangnan.fmme.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.jiangnan.fmme.cell.sim.ReactionProcess() + (...
                    permute(this.rna.molecularWeights(this.rna.matureIndexs(this.aminoacylatedRNAGlobalIndexs))'        * permute(this.freeRNAs,         [1 3 2]),[1 3 2]) + ...
                    permute(this.rna.molecularWeights(this.rna.aminoacylatedIndexs(this.aminoacylatedRNAGlobalIndexs))' * permute(this.aminoacylatedRNAs,[1 3 2]),[1 3 2])) ...
                    / edu.jiangnan.fmme.util.ConstantUtil.nAvogadro;
            end
        end
    end
end

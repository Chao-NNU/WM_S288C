%Transcription
%
% @wholeCellModelID Process_Transcription
% @name             Transcription
% @description
%   Biology
%   ===============
%   Transcription is the first step in the synthesis of functional gene
%   products where RNA polymerase and several accessory enzymes translate
%   transcription units, or regions of the DNA containing 1 or more genes, into
%   RNA molecules. Following transcription the RNA molecules follow one of
%   two pathways:
%   - mRNAs are used as templates for translation (see translation process)
%   - r/s/tRNAs which are transcribed as molecules containing multiple genes are
%     cleavaged into their individual genes (see RNA processing process), are
%     modified at several bases (see RNA modification process) to improve their
%     stability and enhance their catalytic activity, and finally act as
%     ribozymes (r/sRNAs) or as adaptors between mRNAs and the amino acids they
%     code for (tRNAs).
%
%   Transcription begins with the recruitment of RNA polymerase to a promoter
%   with the help of the sigma initiation factor and possiblity transcription
%   factors. Next elongation factors are recruited, RNA begins to be
%   polymerized, and sigma factor is released. Finally the RNA polymerase
%   reaches  a terminator at the the end of the transcription unit, and
%   with the help of termination factors releases the polymerized RNA and
%   dissociates from the DNA. Termination in E. coli occurs via either a
%   Rho-dependent (50%) or Rho-independent mechanism (50%). Rho-dependent
%   termination is catalyzed by the hexameric ATP-dependent helicase Rho.
%   Rho is not essential in B. subtilis [PUB_0234]. Rho-independent
%   termination occur via the intrinsic properties RNA which disrupt RNA
%   polymerase-DNA binding. Terminator hairpins are not predicted in M.
%   genitalium. In E. coli termination is incompetition with
%   antitermination. However antitermination has not been reported in any
%   mcyoplasma [PUB_0182].
%
%   As soon as the RNA begins to polymerized, even prior to termination, the
%   mRNA transcripts may be bound by ribosomes and polymerized. For simplicity,
%   our model doesn't represent this phenomenon.
%
%   Knowledge Base
%   ===============
%   The transcription unit structure was compiled from several sources:
%   - Primary reports of cotranscribed genes [PUB_0176, PUB_0182, PUB_0186,
%     PUB_0188, PUB_0188, PUB_0244, PUB_0247, PUB_0248, PUB_0249, PUB_0251]
%   - OperonDB database of cotranscribed genes [PUB_0250]
%   - Conservation of gene order across multiple species
%   - Related function of adjacent genes
%   - Expression levels as measured by microarrays of adjacent genes [PUB_0569]
%   - Strandedness of adjacent genes
%   - Weiner, Hermann, and Browning computational model of promoters and
%     transcription unit start sites [PUB_0411]
%
%   The transcription unit structure is organized in the knowledge base, and is
%   loaded into this process by the initializeConstants method.
%
%   The knowledge base also contains the measured expression and half-lives of
%   many transcripts. These values are loaded by initializeConstants and fit by
%   simulation.fitConstants to be consistent with other processes.
%
%   Representation
%   ===============
%   substrates, enzymes, boundEnzymes, and RNAs represent the counts of
%   metabolites, free transcription enzymes, transcription enzymes bound to RNAs
%   and RNA polymerases, and nascent RNAs.
%
%   rnaPolymerases.states represents the current state / pseudostate (actively
%   transcribing, specifically bound, non-specifically bound, free,
%   non-existent) of each RNA polymerase, where each state is indicated by the
%   enumeration:
%   - rnaPolymeraseActivelyTranscribingValue
%   - rnaPolymeraseSpecificallyBoundValue
%   - rnaPolymeraseNonSpecificallyBoundValue
%   - rnaPolymeraseFreeValue
%   - rnaPolymeraseNotExistValue (state exists as a way to account for memory
%     allocated for future RNA polymerases)
%   For actively transcribing polymerases rnaPolymerases.states also represents
%   the position of the polymerase along the transcription unit.
%
%   That is entries of rnaPolymerases.states with the following values
%   corresponding to these states:
%     >= RNAPolymerases.activelyTranscribingValue: RNA polymerases position on genome actively transcribing
%     == RNAPolymerases.specificallyBoundValue:    RNA polymerase specifically bound
%     == RNAPolymerases.nonSpecificallyBoundValue: RNA polymerase non-specifically bound
%     == RNAPolymerases.freeValue:                 RNA polymerase free
%     == RNAPolymerases.notExistValue:             RNA polymerase doesn't exist
%
%   transcripts.boundTranscriptionUnits represents the particular transcription
%   unit to which each actively transcribing and specifically bound polymerase
%   is bound.
%
%   rnaPolymeraseStateExpectations represents the expected occupancies of the
%   RNA polymerase states.
%
%   transcriptionUnitBindingProbabilities represents the relative affinity of
%   RNA polymerases for the promoters of each transcription unit.
%   transcriptionFactorBindingProbFoldChange represents the fold change affect of
%   transcription factors on the relative affinities of the RNA polymerse for
%   the promoters. RNA polymerases are assigned to
%   transcription units weighted by the product of
%   transcriptionUnitBindingProbabilities and
%   transcriptionFactorBindingProbFoldChange.
%
%   Initialization
%   ===============
%   All RNAs are initialized to the mature state. This is implemented by the
%   simulation class initializeState method.
%
%   RNA polymerases are initialized to their steady state:
%   - Each RNA polymerase is randomly assigned (with replacement) to one of the
%     actively transcribing, specifically bound, non-specifically bound, or free
%     states weighted by the expected occupancy of each state
%     (rnaPolymeraseStateExpectations)
%   - Actively transcribing and specifically bound polymerases randomly assigned
%     to transcription units weighted by their transcription rates
%     (transcriptionUnitBindingProbabilities)
%   - Each transcription unit to which an actively transcribing polymerase has
%     been assigned is divided into 1 segment for each polymerase
%   - Actively transcribing polymerases randomly assigned to positions within
%     the assigned segment of their assigned transcription unit (positions near
%     the segment border are not allowed to prevent polymerases from being too
%     close to each other) with uniform probably.
%
%   Simulation
%   ===============
%   Evolves the state of RNA polymerase using a markov chain model with four
%   states:
%   - actively translating
%   - specifically bound
%   - non-specifically bound
%   - free
%
%   Transition probabilities are designed to maintain the occupancy of each
%   state within a narrow window around their expected values. Transition
%   probability are determined by four logistic control functions. These can
%   be tuned with the constants
%   - rnaPolymeraseStateExpectations
%
%   RNA polymerase are created in the free state.
%
%   Actively transcribing state:
%   1. Release sigma factor if after first second of elongation
%   2. Elongate transcript according to nucleic acid limits (substrates)
%      if elongation factors are available
%   3. If transcription complete and termination factor available
%      - release transcript
%      - transition RNA polymerase to free state
%      - increment gene expression
%      Otherwise remain in active state
%
%   Specifically bound State:
%   - Can transition to active, specifically bound, non-specifically bound,
%     or free states
%   - Transition into state only if a free sigma factor is available
%     1. Decrement number of free sigma factors
%     2. Pick a transcription unit (tu) to bind to according to
%     Expression transcription unit i~prob(ribosome releases tu i|ribosome active)
%                      =prob(ribosome within RNA polymerase elongation rate bases of length of tu i|ribosome active)
%                      =prob(ribosome within RNA polymerase elongation rate bases of length of tu i|ribosome active, bound to tu i)*prob(ribosome bound to tu i|ribosome active)
%                      =[(RNA polymerase transcription rate)/(length of tu i)] * [(length of tu i)*prob(binding tu i|binding)]
%     prob(binding tu i | binding)~expression tu i
%
%   Non-specifically bound state:
%   - Can transition to specifically bound, non-specifically bound, or free
%     states
%
%   Free state:
%   - Can transition to specifically bound, non-specifically bound, or free
%     states
%
%   Algorithm
%   +++++++++++++++
%   1. Randomly transition RNA polymerases among activlely transcribing,
%      specifically bound, non-specifically, bound, and free states weighted by
%      state transition probabilities. Update rnaPolymerases.states.
%   2. Randomly assign RNA polymerases entering the specifically bound to
%      specific transcription units weighted by the product of
%      transcriptionUnitBindingProbabilities and
%      transcriptionFactorBindingProbFoldChange. Update
%      transcripts.boundTranscriptionUnits.
%   3. Assign RNA polymerase entering the actively transcribing state sigma
%      factors. Update enzymes and boundEnzymes.
%   4. Simulate RNA polymerization by actively transcribing RNA polymerases with
%      the aid of elongation factors. Allocate available nucleic acids among the
%      actively transcribing RNA polymerases. Release sigma factors from RNA
%      polymerases that started at the beginning of the transcription unit and
%      progressed. Update rnaPolymerases.states. Update substrates. Update enzymes
%      and boundEnzymes.
%   5. If termination factors are available dissociate RNA polymerases which
%      have reached the terminus of the transcription they're bound to, and
%      release RNAs. Update rnaPolymerases.states and
%      transcripts.boundTranscriptionUnits. Increment RNAs.
%
%   References
%   ===========
%   1. McClure, W. R. 1985. Mechanism and control of transcription
%      initiation in prokaryotes. Annu. Rev. Biochem. 54:171-204.
%      [PUB_0775]
%   2. Ciampi MS (2006). Rho-dependent terminators and transcription
%      termination. Microbiology. 152(9):2515-28 [PUB_0233].
%   3. Nudler E, Gottesman ME (2002). Transcription termination and
%      anti-termination in E. coli. Genes Cells. 7(8):755-68. [PUB_0662]
%   4. Washio T, Sasayama J, Tomita M (1998). Analysis of complete genomes
%      suggests that many prokaryotes do not rely on hairpin formation in
%      transcription termination. Nucleic Acids Res. 26(23):5456-63
%      [PUB_0234]
%   5. Peterson JD, Umayam LA, Dickinson T, Hickey EK, White O (2001). The
%      Comprehensive Microbial Resource. Nucleic Acids Res. 29(1):123-5.
%      [PUB_0182]
%   6. Shepherd N, Dennis P, Bremer H (2001). Cytoplasmic RNA Polymerase in
%      Escherichia coli. J Bacteriol. 183(8): 2527-34. [PUB_0784]
%   7. Klumpp S, Hwa T (2008). Growth-rate-dependent partitioning of RNA
%      polymerases in bacteria. Proc Natl Acad Sci U S A. 105(21):
%      20245-50. [PUB_0785]
%   8. Grigorova IL, Phleger NJ, Mutalik VK, Gross CA (2006). Insights into
%      transcriptional regulation and sigma competition from an equilibrium
%      model of RNA polymerase binding to DNA. Proc Natl Acad Sci U S A.
%      103(14): 5332-7. [PUB_0786]
%
% Author: Markus Covert, mcovert@stanford.edu
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 1/4/2010
%
%TODO: require Mg2+ or Mn2+ as cofactor for transcription
classdef Transcription < edu.jiangnan.fmme.cell.sim.Process & edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'rnaPolymeraseElongationRate';
            'transcriptionUnitBaseCounts';
            'transcriptionUnitBindingProbabilities';
            };
        fittedConstantNames__      = {   %names of fitted constant properties
            'transcriptionUnitBindingProbabilities'
            };
        localStateNames__          = {   %names of simulation state properties redundant with timecourses in this or other processes or the simulation
            'RNAs'};
    end
    
    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {}; %whole cell model IDs of stimuli
        
        substrateWholeCellModelIDs = { %whole cell model IDs of substrates
            'ATP'; 'CTP'; 'GTP'; 'UTP';
            'AMP'; 'CMP'; 'GMP'; 'UMP';
            'ADP'; 'PPI'; 'H2O'; 'H'};
        substrateIndexs_ntp         = (1:4)'; %indices within substrates of NTPs
        substrateIndexs_nmp         = (5:8)'; %indices within substrates of NMPs
        substrateIndexs_atp         = 1;      %index within substrates of ATP
        substrateIndexs_adp         = 9;      %index within substrates of ADP
        substrateIndexs_diphosphate = 10;     %index within substrates of diphosphate
        substrateIndexs_water       = 11;     %index within substrates of water
        substrateIndexs_hydrogen    = 12;     %index within substrates of hydrogen
		
		enzymeWholeCellModelIDs = {       						%enzyme whole cell model ids
			'YEL009C_MONOMER';			  						%General control protein GCN4
			'YER164W_YGR116W_YNL230C_YPL046C_YPR133C_PENTAMER'; %Transcription elongation factor complex
			'YGL044C_MONOMER';									%Transcription termination factor RNA15
			'YER148W_MONOMER';			  						%TATA-binding protein (TBP)
			'RNA_POLYMERASE_I';									%DNA-directed RNA polymerase I complex
			'RNA_POLYMERASE_II'};								%DNA-directed RNA polymerase II complex
			
		
        enzymeIndexs_transcriptionFactors    = (1:4)'; %indices within enzymes of transcription factors
        enzymeIndexs_GCN4             		 = 1;      %index within enzymes of General control protein GCN4
		enzymeIndexs_elongationFactor        = 2;      %index within enzymes of elongation factor
		enzymeIndexs_terminationFactor       = 3;      %index within enzymes of termination factor RNA15
        enzymeIndexs_antiterminationFactor   = 4;      %index within enzymes of TATA-binding protein (TBP)
		enzymeIndexs_rnaPolymerase1          = 5;      %index within enzymes of DNA-directed RNA polymerase II
		enzymeIndexs_rnaPolymerase2          = 6;      %index within enzymes of DNA-directed RNA polymerase III
		
        complexIndexs_DnaA_ATP                         %indices within protein complexes of DnaA-ATP
        transcriptionUnitIndexs_DnaAR12345Boxes        %indices of transcription units containing the functional DnaA boxes R1-5
		
		TUIndexs									   %indices of transcription units
		
    end
    
    %fixed biological constants
    properties
        rnaPolymeraseElongationRate                   %RNA polymerase elongation rate (50 nucleotides per second per RNA polymerase) [PUB_0562, PUB_0563]
        transcriptionUnitBaseCounts                   %transcription unit base counts
        transcriptionUnitBindingProbabilities	      %transcription unit binding probabilities
        stateTransitionProbabilities                  %transition probabilities among RNA polymerase states
    end
    
    %global state, stored locally for convenience
    properties
        RNAs                                          %copy number of transcription units
    end
    
    %global state (referenced locally for convenience)
    properties
        rnaPolymerases   %RNA Polymerase state class
        transcripts      %New Transcripts state class
    end
    
    %constructor
    methods
        function this = Transcription(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.Process(wholeCellModelID, name);
        end
    end
    
    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.jiangnan.fmme.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.rnaPolymerases = simulation.state('RNAPolymerase');
            this.transcripts = simulation.state('Transcript');
            
            this.states = [this.states; {this.rnaPolymerases; this.transcripts}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            %call super class method
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.Process(knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect(knowledgeBase, simulation, varargin{:});
			for i = 1:16
				g{i} = knowledgeBase.genomes(i);

				this.transcriptionUnitBaseCounts{i} = reshape([g{i}.transcriptionUnits.baseCount],[],length(g{i}.transcriptionUnits))';
				this.transcriptionUnitBaseCounts{i} = this.transcriptionUnitBaseCounts{i}(:, this.substrateMetaboliteGlobalIndexs);
				%this.transcriptionUnitBaseCounts = reshape([knowledgeBase.transcriptionUnits.baseCount],[],length(knowledgeBase.transcriptionUnits))';
				%this.transcriptionUnitBaseCounts = this.transcriptionUnitBaseCounts(:, this.substrateMetaboliteGlobalIndexs);
				for j = 1:length(g{i}.transcriptionUnits)
					this.transcriptionUnitBindingProbabilities{i}(j) = zeros(size(g{i}.transcriptionUnits(j).sequenceLength));
					%this.transcriptionUnitBindingProbabilities{i}  = zeros(size(this.transcripts.transcriptionUnitLengths));
				end
			end
			this.complexIndexs_DnaA_ATP = find(ismember({knowledgeBase.proteinComplexs.wholeCellModelID}', {
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
                'YML065W_YBR060C_YLL004W_YPR162C_YNL261W_YHR118C_6MER_6ATP'})); %Orc-ATP 6mer	
			for i = 1:16
			[~, idxs{i}] = ismember({
				'ARS consensus sequence 1'
				'ARS consensus sequence 2'
				'ARS consensus sequence 3'
				'ARS consensus sequence 4'
				'ARS consensus sequence 5'
				}, {knowledgeBase.genomes(i).genomeFeatures.name});
				
			tus{i} = [knowledgeBase.genomes(i).genomeFeatures(idxs{i}).transcriptionUnits];
			tu{i}  = min([tus{i}.idx])';
			this.transcriptionUnitIndexs_DnaAR12345Boxes{i} = find([knowledgeBase.genomes(i).transcriptionUnits.idx]' ==  tu{i});
			end
			
			for i = 1:6447
				B(i) = knowledgeBase.genes(i).genomes.idx;
			end
			for i = 1:16
				this.TUIndexs{i} = find(B==i);
			end
        end
        
        %retrieve state from simulation
        function copyFromState(this)
            this.copyFromState@edu.jiangnan.fmme.cell.sim.Process();
            
            this.RNAs = this.rna.counts(this.rna.nascentIndexs, this.compartment.cytosolIndexs, :);
        end
        
        %send state to simulation
        function copyToState(this)
            this.copyToState@edu.jiangnan.fmme.cell.sim.Process();
            this.rna.counts(this.rna.nascentIndexs, this.compartment.cytosolIndexs, :) = this.RNAs;
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemoryForState(this, numTimePoints)
            this.allocateMemoryForState@edu.jiangnan.fmme.cell.sim.Process(numTimePoints);
            numTranscriptionUnits = 6447;%manual added
            this.RNAs = zeros(numTranscriptionUnits, 1, numTimePoints);
        end
    end
    
    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, states)
            %% initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %% substrate and byproducts
            %Contributions, by process, to RNA metabolism:
            %                 NTPs + H2O --(transcription)--> nascent RNA + PPi + H
            %          nascent RNA + H2O --(processing)-----> processed RNA + NMPs
            %processed RNA + metabolites --(modification)---> modified RNA + metabolites
            %               modified RNA --(decay)----------> NMPs + modified NMPs
            %               NMPs + 2 ATP --(charging)-------> NTPs + 2 ADP
            %        modified NMPs + H2O --(catabolism)-----> modified bases + Pi
            %
            %Here we only account for the contributions of transcription.
            invMat = edu.jiangnan.fmme.util.ComputationUtil.invertCompositionMatrix(this.rna.nascentRNAMatureRNAComposition);
			
			transcribedNTPss = zeros(4,6447);
			
			for i = 1:16
			transcribedNTP{i} = this.transcriptionUnitBaseCounts{i}(:, this.substrateIndexs_nmp)';
			transcribedNTPss(:,this.TUIndexs{i})     = transcribedNTP{i};
			end
			
			transcribedNTPs = transcribedNTPss * invMat * states.rnaProductions;
			
            %transcribedNTPs = this.transcriptionUnitBaseCounts(:, this.substrateIndexs_nmp)' * invMat * states.rnaProductions;
            bmProd(this.substrateIndexs_ntp)         = bmProd(this.substrateIndexs_ntp)         + transcribedNTPs;
            byProd(this.substrateIndexs_diphosphate) = byProd(this.substrateIndexs_diphosphate) + sum(transcribedNTPs);
            bmProd(this.substrateIndexs_water)       = bmProd(this.substrateIndexs_water)       + sum(states.rnaProductions);
            byProd(this.substrateIndexs_hydrogen)    = byProd(this.substrateIndexs_hydrogen)    + sum(states.rnaProductions);
            
            %% enzymes
            %RNA polymerase
			len = this.transcripts.transcriptionUnitLengths;
			lens = zeros(6447,1);
			for i = 1:16
			lens(this.TUIndexs{i},1) = len{i};
			end
			
            fractionInitiating = this.rnaPolymeraseElongationRate * sum(invMat * states.rnaProductions0) / ... %fraction of active RNA polymerases that are initiating
                (lens' * invMat * states.rnaProductions0);
            minEnzExp(this.enzymeIndexs_rnaPolymerase1) = 4.5867 * ...
                (lens' * invMat * states.rnaProductions0) / ...
                this.rnaPolymeraseElongationRate / ...
                this.rnaPolymerases.stateExpectations(this.rnaPolymerases.activelyTranscribingIndex) / ...
                (1 - fractionInitiating);
            
            %sigma, elongation, termination factors
            minEnzExp(this.enzymeIndexs_GCN4) = minEnzExp(this.enzymeIndexs_rnaPolymerase1) * ...
                (this.rnaPolymerases.stateExpectations(this.rnaPolymerases.activelyTranscribingIndex) * fractionInitiating + ...
                this.rnaPolymerases.stateExpectations(this.rnaPolymerases.specificallyBoundIndex));
            minEnzExp(this.enzymeIndexs_elongationFactor)  = 2;
            minEnzExp(this.enzymeIndexs_terminationFactor) = 2;
        end
        
        %Initialize RNA polymerase, bound transcription units, transcription factor states
        %Assumptions:
        %- Metabolic cost of the transcription of these nucleic acids is negligible
        %- Probability that so many RNA polymerases bind a given transcription unit
        %  that they can't pack along the transcription unit without steric
        %  interactions is negligible. Thus we don't try to handle this case
        %  separately
        function initializeState(this)
            rnaPols = this.rnaPolymerases;
            trnspts = this.transcripts;

            %% states
			for i = 1:16
            trnsptStrnds{i} = 2 - trnspts.transcriptionUnitDirections{i};
            trnsptLens{i} = trnspts.transcriptionUnitLengths{i};
            trnsptDirs{i} = 2 * trnspts.transcriptionUnitDirections{i} - 1;
            trnspt5Coords{i} = trnspts.transcriptionUnitFivePrimeCoordinates{i};
            trnsptStarts{i} = ...
                trnspt5Coords{i} - (1 - trnspts.transcriptionUnitDirections{i}) .* (trnspts.transcriptionUnitLengths{i}-1);
				
            trnsptStarts{i}(trnsptDirs{i} == 1) = ...
                trnsptStarts{i}(trnsptDirs{i} == 1) - ...
                this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase1);
            trnsptStarts{i}(trnsptDirs{i} == -1) = ...
                trnsptStarts{i}(trnsptDirs{i} == -1) + ...
                this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase1);
			end
			
			%% enzymes
            freTFs = this.enzymes(this.enzymeIndexs_transcriptionFactors);
            bndTFs = this.boundEnzymes(this.enzymeIndexs_transcriptionFactors);
            frePls = this.enzymes(this.enzymeIndexs_rnaPolymerase1);
            bndPls = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase1);
            freHls = this.enzymes(this.enzymeIndexs_rnaPolymerase2);
            bndHls = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase2);
            pls = frePls + bndPls + freHls + bndHls;
			
			%keyboard
			for i = 1:16
            %% allocate
			pl(i) = fix(pls * (length(trnspts.transcriptionUnitWholeCellModelIDs{i}) / 6401));
            rnaPols.states{i} = repmat(rnaPols.notExistValue, [2 * pl(i), 1, 1]);
            rnaPols.positionStrands{i} = zeros([2 * pl(i), 2, 1]);
            trnspts.boundTranscriptionUnits{i} = repmat(trnspts.nullTranscriptValue, [2 * pl(i), 1, 1]);
            trnspts.boundTranscriptProgress{i} = repmat(trnspts.nullTranscriptValue, [2 * pl(i), 1, 1]);
            trnspts.boundTranscriptChromosome{i} = repmat(trnspts.nullTranscriptValue, [2 * pl(i), 1, 1]);			
            end
			%keyboard
            %% initialize
			for i = 1:16
				pState    = rnaPols.stateExpectations;
				probs{i}  = this.computeRNAPolymeraseTUBindingProbabilities(i);			
				pBinds{i} = trnspts.transcriptionUnitLengths{i} .* probs{i}(:, 1);
				pSpBinds{i} = pBinds{i};
			end
			for i = 16
            %for j = 1:pls %change from i to j
			for j = 1:pl(i)
				if frePls == 0 || ~any(pBinds{i}) || freTFs(this.enzymeIndexs_elongationFactor) == 0
                    pState(rnaPols.activelyTranscribingIndex) = 0;
                end				
				if frePls == 0 || ~any(pSpBinds{i}) || freTFs(this.enzymeIndexs_GCN4) == 0
                    pState(rnaPols.specificallyBoundIndex) = 0;
                end				
				if frePls == 0
                    pState(rnaPols.nonSpecificallyBoundIndex) = 0;
                end				
                if ~any(pState)
                    pState(rnaPols.freeIndex) = 1;
                end
				
                %randomly select state
				state{i} = this.randStream.randsample(4, 1, true, pState);

                %Actively transcribing, specifically bound: randomly select
                %transcription unit
                switch state{i}	
                    case rnaPols.activelyTranscribingIndex	%1
                        tmp_pBinds{i} = pBinds{i};
                        posStrnds{i} = zeros(0, 2);
						while any(tmp_pBinds{i})
                            iTU{i} = this.randStream.randsample(numel(pBinds{i}), 1, true, tmp_pBinds{i});
                            tmp_pBinds{i}(iTU{i}) = 0;
							%keyboard
                            [~, posStrnds{i}] = this.bindProteinToChromosomeStochastically(i, ...
                                this.enzymeIndexs_rnaPolymerase1, 1, [trnsptStarts{i}(iTU{i})+1 trnsptStrnds{i}(iTU{i})], trnsptLens{i}(iTU{i})-2);
						if ~isempty(posStrnds{i})
                                break;
                            end
                        end
						%keyboard
						if isempty(posStrnds{i})
                            %throw(MException('Transcription:error', 'unable to bind chromosome'));
							warning('unable to bind chromosome');
                        end
						
                        posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) = posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase1);
                        posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) = posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase1);
						
                        if trnspts.transcriptionUnitDirections{i}(iTU{i})
                            state{i} = posStrnds{i}(:, 1) - trnspt5Coords{i}(iTU{i}) + 1;
                        else
                            state{i} = trnspt5Coords{i}(iTU{i}) - posStrnds{i}(:, 1) + 1;
                        end
						
						%keyboard
                        rnaPols.states{i}(j) = state{i};
                        rnaPols.positionStrands{i}(j, :) = posStrnds{i};
                        trnspts.boundTranscriptProgress{i}(j) = state{i};
                        trnspts.boundTranscriptionUnits{i}(j) = iTU{i};
                        trnspts.boundTranscriptChromosome{i}(j) = 1; %initialize to first chromosome
                        
                        frePls = frePls - 1;
                        bndPls = bndPls + 1;
						
                    case rnaPols.specificallyBoundIndex	%2
                        while any(pSpBinds{i})

                            iTU{i} = this.randStream.randsample(numel(pSpBinds{i}), 1, true, pSpBinds{i});
                            pSpBinds{i}(iTU{i}) = 0;

                            [~, ~, posStrnds{i}] = this.bindProteinToChromosome(i,[trnspt5Coords{i}(iTU{i})-trnsptDirs{i}(iTU{i}) trnsptStrnds{i}(iTU{i})], this.enzymeIndexs_rnaPolymerase2, 1);			
                            if ~isempty(posStrnds{i})
                                break;
                            end
                        end
						
                        if isempty(posStrnds{i})
                            %throw(MException('Transcription:error', 'unable to bind chromosome'));
							warning('unable to bind chromosome');
                        end
                        
                        rnaPols.states{i}(j) = rnaPols.specificallyBoundValue;
                        rnaPols.positionStrands{i}(j, :) = posStrnds{i};
                        trnspts.boundTranscriptProgress{i}(j) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptionUnits{i}(j) = iTU{i};
                        trnspts.boundTranscriptChromosome{i}(j) = 1; %initialize to first chromosome
                        
                        bndHls = bndHls + 1;
                        frePls = frePls - 1;
                        freTFs(this.enzymeIndexs_GCN4) = freTFs(this.enzymeIndexs_GCN4) - 1;
						
                    case rnaPols.nonSpecificallyBoundIndex	%3
					
                        [~, posStrnds{i}] = this.bindProteinToChromosomeStochastically(i, this.enzymeIndexs_rnaPolymerase1, 1);
						
                        if isempty(posStrnds{i})
                            %throw(MException('Transcription:error', 'unable to bind chromosome'));
							warning('unable to bind chromosome');
                        end
						
                        posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) = posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase1);
                        posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) = posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase1);
                        
                        rnaPols.states{i}(j) = rnaPols.nonSpecificallyBoundValue;
						
                        rnaPols.positionStrands{i}(j, :) = posStrnds{i};
                        trnspts.boundTranscriptionUnits{i}(j) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptProgress{i}(j) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptChromosome{i}(j) = trnspts.nullTranscriptValue;
                        
                        frePls = frePls - 1;
                        bndPls = bndPls + 1;						
                    case rnaPols.freeIndex	%4
                        rnaPols.states{i}(j) = rnaPols.freeValue;
                        rnaPols.positionStrands{i}(j, :) = 0;
                        trnspts.boundTranscriptionUnits{i}(j) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptProgress{i}(j) = trnspts.nullTranscriptValue;
                        trnspts.boundTranscriptChromosome{i}(j) = trnspts.nullTranscriptValue;
					end	
                end
            end
            
            %% store states of enzymes
            this.enzymes(this.enzymeIndexs_transcriptionFactors)      = freTFs;
            this.boundEnzymes(this.enzymeIndexs_transcriptionFactors) = bndTFs;
            this.enzymes(this.enzymeIndexs_rnaPolymerase1)            = frePls;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase1)       = bndPls;
            this.enzymes(this.enzymeIndexs_rnaPolymerase2)            = freHls;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase2)       = bndHls;
            
            %% decrement counts of RNAs
            this.copyToState();
            
            comp = this.compartment;
            rna = this.rna;
            
            matureRnaWt = max(0, sum(rna.dryWeight) - trnspts.dryWeight);
            initRnaCnts = rna.counts;
            while sum(rna.dryWeight) > matureRnaWt
                idx = this.randStream.randsample(size(rna.counts, 1), 1, false, rna.counts(:, comp.cytosolIndexs));
                rna.counts(idx, comp.cytosolIndexs) = rna.counts(idx, comp.cytosolIndexs) - 1;
            end
            rna.counts(:, comp.cytosolIndexs) = ...
                + rna.counts(:, comp.cytosolIndexs) ...
                + rna.updateExternalState(-(initRnaCnts(:, comp.cytosolIndexs) - rna.counts(:, comp.cytosolIndexs)), true);
            
            this.copyFromState();
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
            result(this.substrateIndexs_ntp) = 2 * ...
                sum(this.metabolite.nmpComposition, 2) * sum(cell2mat(this.rnaPolymerases.nActive))* ...
                this.rnaPolymeraseElongationRate;
            result(this.substrateIndexs_water) = sum(cell2mat(this.rnaPolymerases.nActive));
        end
        
        %simulation
        function evolveState(this)
            %% define and allocate variables
            
            %states
            c = this.chromosome;
            rnaPols = this.rnaPolymerases;
            trnspts = this.transcripts;
			
            %numbers of enzymes
            nFreTfs  = this.enzymes(this.enzymeIndexs_transcriptionFactors);
            nBndTFs  = this.boundEnzymes(this.enzymeIndexs_transcriptionFactors);
            nFrePols = this.enzymes(this.enzymeIndexs_rnaPolymerase1);
            nBndPols = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase1);
            nFreHols = this.enzymes(this.enzymeIndexs_rnaPolymerase2);
            nBndHols = this.boundEnzymes(this.enzymeIndexs_rnaPolymerase2);
            nTotPols = nFrePols + nBndPols + nFreHols + nBndHols;

			for i = 16			
            %properties of RNA polymerases
            pStTrns{i} = this.stateTransitionProbabilities{i}; %probabilities of state transitions
            
            %transcription unit properties
            tuDirs{i} = 2 * trnspts.transcriptionUnitDirections{i} - 1;
            tuStrnds{i} = c.transcriptionUnitStrands{i};
            tuLens{i} = trnspts.transcriptionUnitLengths{i};
            tu5Coords{i} = trnspts.transcriptionUnitFivePrimeCoordinates{i};
            tuSeqs{i} = trnspts.transcriptionUnitSequences{i};  %sequences
            
            %% Update states of RNA polymerases
            %allocate space to store states of new RNA polymerases
			if nTotPols > size(rnaPols.states{i}, 1)
                rnaPols.states{i} = [
                    rnaPols.states{i};
                    rnaPols.notExistValue(ones(2 * nTotPols - size(rnaPols.states{i}, 1), 1), 1)];
                rnaPols.positionStrands{i} = [
                    rnaPols.positionStrands{i};
                    zeros(2 * nTotPols - size(rnaPols.positionStrands{i}, 1), 2)];
                trnspts.boundTranscriptionUnits{i} = [
                    trnspts.boundTranscriptionUnits{i};
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome{i}, 1), 1), 1)];
                trnspts.boundTranscriptProgress{i} = [
                    trnspts.boundTranscriptProgress{i};
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome{i}, 1), 1), 1)];
                trnspts.boundTranscriptChromosome{i} = [
                    trnspts.boundTranscriptChromosome{i};
                    trnspts.nullTranscriptValue(ones(2 * nTotPols - size(trnspts.boundTranscriptChromosome{i}, 1), 1), 1)];
            end			
			
            %indices of RNA polymerases in various states
            actPols{i} = find(rnaPols.states{i} >= rnaPols.activelyTranscribingValue);
            sbPols{i} = find(rnaPols.states{i} == rnaPols.specificallyBoundValue);
            freeNsbPols{i} = find(...
                rnaPols.states{i} == rnaPols.nonSpecificallyBoundValue | ...
                rnaPols.states{i} == rnaPols.freeValue);
            sbPols{i} = sbPols{i}(this.randStream.randperm(numel(sbPols{i})));
            freeNsbPols{i} = freeNsbPols{i}(this.randStream.randperm(numel(freeNsbPols{i})));
            
            %number of new RNA polymerases
			
            nNewPols{i} = nTotPols - (numel(actPols{i}) + numel(sbPols{i}) + numel(freeNsbPols{i}));			
            
            %dissociate any free RNA polymerase holoenzymes (created say by
            %their release from chromosome by another protein) into RNA
            %polymerase core and sigma factor
			nFrePols = nFrePols + nFreHols;
            nFreTfs(this.enzymeIndexs_GCN4) = nFreTfs(this.enzymeIndexs_GCN4) + nFreHols;
            nFreHols = 0;			
            %initiating specifically bound polymerases
			
            nInitMax{i} = this.randStream.stochasticRound(numel(sbPols{i}) * pStTrns{i}(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex));
			%#ok<*PROP>
            if nInitMax{i} > 0 %meet the condition
                releasedProteins{i} = this.releaseProteinFromSites(i, rnaPols.positionStrands{i}(sbPols{i}(1:nInitMax{i}), :), false, this.enzymeIndexs_rnaPolymerase2, true, true);
				if ~isequal(releasedProteins{i}(this.enzymeIndexs_rnaPolymerase2), nInitMax{i})
                    %throw(MException('Transcription:error', 'Unable to release protein'));
					warning('Unable to release protein');
                end
				
                %try to initiate RNA polymerases
                positions{i} = tu5Coords{i}(trnspts.boundTranscriptionUnits{i}(sbPols{i}(1:nInitMax{i})));
                tfs{i} = this.bindProteinToChromosome(i, [positions{i} rnaPols.positionStrands{i}(sbPols{i}(1:nInitMax{i}), 2)], ...
                    this.enzymeIndexs_rnaPolymerase2);
                rnaPols.states{i}(sbPols{i}(tfs{i})) = rnaPols.activelyTranscribingValue;
                trnspts.boundTranscriptProgress{i}(sbPols{i}(tfs{i})) = trnspts.newTranscriptValue;
                rnaPols.positionStrands{i}(sbPols{i}(tfs{i}), 1) = positions{i}(tfs{i}, 1);
                
                %rebind RNA polymerases that couldn't move forward
                
				if ~all(this.bindProteinToChromosome(i, rnaPols.positionStrands{i}(sbPols{i}(~tfs{i}), :), ...
                        this.enzymeIndexs_rnaPolymerase2))
                    %throw(MException('Transcription:error', 'Unable to unbind protein'));
					warning('Unable to unbind protein')
                end
            end
			
            %specifically bound polymerases becoming non-specifically bound / free
            nUnsb{i} = max(0, numel(sbPols{i}) - nInitMax{i} - this.randStream.stochasticRound(numel(sbPols{i}) * ...
                pStTrns{i}(rnaPols.specificallyBoundIndex, rnaPols.specificallyBoundIndex)));
            if nUnsb{i} > 0 %not meet the condition
                idxs{i} = sbPols{i}(end-nUnsb{i}+1:end);
                posStrnds{i} = rnaPols.positionStrands{i}(idxs{i}, :);
                releasedProteins{i} = this.releaseProteinFromSites(i, posStrnds{i}, false, this.enzymeIndexs_rnaPolymerase2, true, true);
                
				if ~isequal(releasedProteins(this.enzymeIndexs_rnaPolymerase2), nUnsb{i})
                    %throw(MException('Transcription:error', 'Unable to unbind protein'));
					warning('Unable to unbind protein');
                end                
				
                rnaPols.states{i}(idxs{i}) = rnaPols.freeValue;
                rnaPols.positionStrands{i}(idxs{i}, :) = 0;
                trnspts.boundTranscriptionUnits{i}(idxs{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress{i}(idxs{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(idxs{i}) = trnspts.nullTranscriptValue;
                nBndHols = nBndHols - nUnsb{i};                
                nFrePols = nFrePols + nUnsb{i};
                nFreTfs(this.enzymeIndexs_GCN4) = nFreTfs(this.enzymeIndexs_GCN4) + nUnsb{i};	
            end
			
            nFreeNsbPols{i} = this.randStream.stochasticRound(min([
                numel(freeNsbPols{i}) * pStTrns{i}(rnaPols.specificallyBoundIndex, rnaPols.freeIndex)
                nFreTfs(this.enzymeIndexs_GCN4)]));		
				
            if nFreeNsbPols{i} > 0 %not meet the condition
                %unbind from non-specifically bound site
                idxs{i} = freeNsbPols{i}(1:nFreeNsbPols{i});
				
				nNsbPols{i} = sum(rnaPols.states{i}(idxs{i}) == rnaPols.nonSpecificallyBoundValue);
                releasedProteins{i} = this.releaseProteinFromSites(i, rnaPols.positionStrands{i}(idxs{i}, :), false, this.enzymeIndexs_rnaPolymerase1, true, true);
				
				if ~isequal(releasedProteins{i}(this.enzymeIndexs_rnaPolymerase1), nNsbPols{i})
                    %throw(MException('Transcription:error', 'Unable to release protein'));
					warning('Unable to release protein');
                end
                
                rnaPols.states{i}(idxs{i}) = rnaPols.freeValue;
                rnaPols.positionStrands{i}(idxs{i}, :) = 0;
                trnspts.boundTranscriptionUnits{i}(idxs{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress{i}(idxs{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(idxs{i}) = trnspts.nullTranscriptValue;
				
                %nFrePols = nFrePols + nNsbPols{i};
                %nBndPols = nBndPols - nNsbPols{i};
                nFrePols = nFrePols + nNsbPols{i};
                nBndPols = nBndPols - nNsbPols{i};				
				
                
                %bind to new specifically bound site
                pBnds{i} = this.computeRNAPolymeraseTUBindingProbabilities(i);
				
                if nnz(c.polymerizedRegions{i}) == 2
                    [~, iTU{i}, posStrnds{i}, nFreeNsbPols{i}] = this.bindProteinToChromosome(i, ...
                        [tu5Coords{i}-tuDirs{i} tuStrnds{i}], ...
                        this.enzymeIndexs_rnaPolymerase2, nFreeNsbPols{i}, pBnds{i}(:, 1), ...
                        true, true, 1, false, []);
                    iChr = 1;
                else
                    [~, iTU{i}, posStrnds{i}, nFreeNsbPols{i}] = this.bindProteinToChromosome(i, ...
                        [tu5Coords{i}-tuDirs{i} tuStrnds{i}
                        tu5Coords{i}-tuDirs{i} tuStrnds{i} + 2], ...
                        this.enzymeIndexs_rnaPolymerase2, nFreeNsbPols{i}, pBnds{i}(:), ...
                        true, true, 1, false, []);
                    iTU{i} = mod(iTU{i} - 1, numel(tuLens{i})) + 1;
                    iChr = ceil(posStrnds{i}(:, 2) / 2);
                end
                
                idxs{i} = freeNsbPols{i}(1:nFreeNsbPols{i});
                rnaPols.states{i}(idxs{i}) = rnaPols.specificallyBoundValue;
                rnaPols.positionStrands{i}(idxs{i}, :) = posStrnds{i};
                trnspts.boundTranscriptionUnits{i}(idxs{i}) = iTU{i};
                trnspts.boundTranscriptProgress{i}(idxs{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(idxs{i}) = iChr;
                
                nBndHols = nBndHols + nFreeNsbPols{i};
                nFrePols = nFrePols - nFreeNsbPols{i};
                nFreTfs(this.enzymeIndexs_GCN4) = nFreTfs(this.enzymeIndexs_GCN4) - nFreeNsbPols{i};
            end
			
            %free polymerases/non-specifically bound becoming free/non-specifically bound
            nsbPols{i} = find(rnaPols.states{i} == rnaPols.nonSpecificallyBoundValue);
			
            if ~isempty(nsbPols{i})
                %set non-specifically bound to free
                releasedProteins{i} = this.releaseProteinFromSites(i, rnaPols.positionStrands{i}(nsbPols{i}, :),...
									  false, this.enzymeIndexs_rnaPolymerase1, true, true);
									  
                if ~isequal(releasedProteins{i}(this.enzymeIndexs_rnaPolymerase1), numel(nsbPols{i}))
                    %throw(MException('Transcription:error', 'Unable to release protein'));
					warning('Unable to release protein');
                end
                rnaPols.states{i}(nsbPols{i}) = rnaPols.freeValue;
                rnaPols.positionStrands{i}(nsbPols{i}, :) = 0;
                trnspts.boundTranscriptionUnits{i}(nsbPols{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress{i}(nsbPols{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(nsbPols{i}) = trnspts.nullTranscriptValue;
                
                %nFrePols = nFrePols + numel(nsbPols{i});
                %nBndPols = nBndPols - numel(nsbPols{i});
                nFrePols = nFrePols + numel(nsbPols{i});
                nBndPols = nBndPols - numel(nsbPols{i});				
            end
			
            freePols{i} = find(rnaPols.states{i} == rnaPols.freeValue);
            %{
			nNsbPols{i} = this.randStream.stochasticRound(min([...
                nFrePols
                numel(freePols{i}) * pStTrns{i}(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) / (pStTrns{i}(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) + pStTrns{i}(rnaPols.freeIndex, rnaPols.freeIndex))]));
			%}
			nNsbPols{i} = this.randStream.stochasticRound(min([...
                nFrePols
                numel(freePols{i}) * pStTrns{i}(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) / (pStTrns{i}(rnaPols.nonSpecificallyBoundIndex, rnaPols.freeIndex) + pStTrns{i}(rnaPols.freeIndex, rnaPols.freeIndex))]));
				
            if nNsbPols{i} > 0 %meet the condition
                %set some to non-specifically bound state
                [nNsbPols{i} posStrnds{i}] = this.bindProteinToChromosomeStochastically(i, this.enzymeIndexs_rnaPolymerase1, nNsbPols{i});
                posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) = posStrnds{i}( isodd(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints5Prime(this.enzymeIndexs_rnaPolymerase1);
                posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) = posStrnds{i}(iseven(posStrnds{i}(:, 2)), 1) + this.enzymeDNAFootprints3Prime(this.enzymeIndexs_rnaPolymerase1);
				                
                rnaPols.states{i}(freePols{i}(1:nNsbPols{i})) = rnaPols.nonSpecificallyBoundValue;
                rnaPols.positionStrands{i}(freePols{i}(1:nNsbPols{i}), :) = posStrnds{i};
                trnspts.boundTranscriptionUnits{i}(freePols{i}(1:nNsbPols{i})) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress{i}(freePols{i}(1:nNsbPols{i})) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(freePols{i}(1:nNsbPols{i})) = trnspts.nullTranscriptValue;
                
                %nFrePols = nFrePols - nNsbPols{i};
                %nBndPols = nBndPols + nNsbPols{i};
				nFrePols = nFrePols - nNsbPols{i};
                nBndPols = nBndPols + nNsbPols{i};
            end
	
            %set newly created RNA polymerases to free state
            if nNewPols{i} > 0
                newPols{i} = find(rnaPols.states{i} == rnaPols.notExistValue, nNewPols{i}, 'first');
                rnaPols.states{i}(newPols{i}) = rnaPols.freeValue;
                rnaPols.positionStrands{i}(newPols{i}, :) = 0;
                trnspts.boundTranscriptionUnits{i}(newPols{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptProgress{i}(newPols{i}) = trnspts.nullTranscriptValue;
                trnspts.boundTranscriptChromosome{i}(newPols{i}) = trnspts.nullTranscriptValue;
            end
		
            %% Transcribe active sequences
            %compute progress of transcription and cost
            usedNTPs = zeros(4, 1);
            if ~isempty(actPols{i}) && nFreTfs(this.enzymeIndexs_elongationFactor) && this.rnaPolymeraseElongationRate > 0
                posStrnds{i} = rnaPols.positionStrands{i}(actPols{i}, :);
                iTUs{i} = trnspts.boundTranscriptionUnits{i}(actPols{i});
                
                %temporarily release RNA polymerase from old position strands
                releasedProteins{i} = ...
                    this.releaseProteinFromSites(i, posStrnds{i}, false, this.enzymeIndexs_rnaPolymerase1, true, true) + ...
                    this.releaseProteinFromSites(i, posStrnds{i}, false, this.enzymeIndexs_rnaPolymerase2, true, true);
                if ~isequal(releasedProteins{i}(this.enzymeIndexs_rnaPolymerase1) + releasedProteins{i}(this.enzymeIndexs_rnaPolymerase2), numel(actPols{i}))
                    %throw(MException('Transcription:error', 'Unable to unbind protein'));
					warning('Unable to unbind protein');
                end
                %extract active sequences (part of sequence that should be transcribed)
                actSeqs{i} = cell(size(actPols{i})); %active sequences (those can be transcribed)
				
                for idx = 1:numel(actPols{i})
                    j = actPols{i}(idx);
                    left{i} = rnaPols.states{i}(j);
                    actSeqs{i}{idx} = tuSeqs{i}{trnspts.boundTranscriptionUnits{i}(j)}(left{i}:end);
                end
                actSeqs{i} = char(actSeqs{i});
                actSeqs{i} = actSeqs{i}(:, 1:min(end, this.rnaPolymeraseElongationRate));
                
                %elongation limits by RNA polymerase kinetic rate and transcription unit extent
                elngMax{i} = min(this.rnaPolymeraseElongationRate, ...
                    tuLens{i}(trnspts.boundTranscriptionUnits{i}(actPols{i})) - rnaPols.states{i}(actPols{i}) + 1);

					%elongation limits by RNA-polymerase / RNA-polymerase interactions
                for strnd = 1 : 4
                    idxs{i} = find(posStrnds{i}(:, 2) == strnd);
                    if numel(idxs{i}) <= 1
                        continue;
                    end
                    [~, order{i}] = sort(posStrnds{i}(idxs{i}, 1));
					
                    if isodd(strnd)
                        elngMax{i}(idxs{i}(order{i})) = min(...
                            elngMax{i}(idxs{i}(order{i})), ...
                            abs(diff([posStrnds{i}(idxs{i}(order{i}), 1); posStrnds{i}(idxs{i}(order{i}(1)), 1) + c.sequenceLen(i)]) - ...
                            this.enzymeDNAFootprints(this.enzymeIndexs_rnaPolymerase1)));
						A = find(elngMax{i} == 0);
						for j = 1:length(A)	
							elngMax{i}(A(j)) = 50;	
						end
                    else
                        elngMax{i}(idxs{i}(order{i})) = min(...
                            elngMax{i}(idxs{i}(order{i})), ...
                            abs(diff([posStrnds{i}(idxs{i}(order{i}(end)), 1) - c.sequenceLen(i); posStrnds{i}(idxs{i}(order{i}), 1)]) - ...
                            this.enzymeDNAFootprints(this.enzymeIndexs_rnaPolymerase1)));
						B = find(elngMax{i} == 0);
						for j = 1:length(B)	
							elngMax{i}(B(j)) = 50;	
						end		
                    end
                end

                %elongation limits by RNA-polymerase DNA interactions (Note: independent of other RNA polymerases)
                [~, ~, ~, tmp_elngMax{i}] = c.isRegionAccessible(i, ...
                    posStrnds{i}, (elngMax{i} + 1) .* tuDirs{i}(iTUs{i}), ...
                    [], this.enzymeGlobalIndexs(this.enzymeIndexs_rnaPolymerase1), true, [], true, false);
				
				if ~all(tmp_elngMax{i})
                    %throw(MException('Transcription:error', 'Unable to bind proteins'));
					warning('Unable to bind proteins');
                end
				
                %elngMax{i} = min(abs(tmp_elngMax{i}) - 1, elngMax{i});
				
				elngMax{i} = abs(min(abs(tmp_elngMax{i}) - 1, elngMax{i}));
				
                %elongation limits by NTPs
                for j = 1:size(actSeqs{i}, 1)
                    actSeqs{i}(j, elngMax{i}(j) + 1:end) = ' ';
                end
                [elngProg, this.substrates(this.substrateIndexs_ntp), usedNTPs] = ...
                    edu.jiangnan.fmme.cell.sim.util.polymerize(...
                    actSeqs{i}, this.substrates(this.substrateIndexs_ntp), 'ACGU', ' ', 0, 0, this.randStream);
                
                %release sigma factor after initiation
                nInits{i} = sum(...
                    rnaPols.states{i}(actPols{i}) == rnaPols.activelyTranscribingValue & ...
                    elngProg > 0);
					
                %nBndHols = nBndHols - nInits{i};
                %nFreTfs(this.enzymeIndexs_GCN4) = nFreTfs(this.enzymeIndexs_GCN4) + nInits{i};
                %nBndPols = nBndPols + nInits{i};
				nBndHols = nBndHols - nInits{i};
                nFreTfs(this.enzymeIndexs_GCN4) = nFreTfs(this.enzymeIndexs_GCN4) + nInits{i};
				nBndPols(i) = nBndPols(i) + nInits{i};
                
                %rebind RNA polymerase to new positions
                rnaPols.states{i}(actPols{i}) = rnaPols.states{i}(actPols{i}) + elngProg;
                rnaPols.positionStrands{i}(actPols{i}, 1) = posStrnds{i}(:, 1) + elngProg .* tuDirs{i}(iTUs{i});
                trnspts.boundTranscriptProgress{i}(actPols{i}) = rnaPols.states{i}(actPols{i});
                
                initPols{i} = rnaPols.states{i}(actPols{i}) == rnaPols.activelyTranscribingValue;
                
				if ...
                        ~all(this.bindProteinToChromosome(i, rnaPols.positionStrands{i}(actPols{i}(~initPols{i}), :), this.enzymeIndexs_rnaPolymerase1)) || ...
                        ~all(this.bindProteinToChromosome(i, rnaPols.positionStrands{i}(actPols{i}( initPols{i}), :), this.enzymeIndexs_rnaPolymerase2))
                    %throw(MException('Transcription:error', 'Unable to bind proteins'));
					warning('Unable to bind proteins');
                end
            end
            %release polymerases which have completed transcription, and increment
            %expression of the corresponding transcription unit
            if ~isempty(actPols{i}) && nFreTfs(this.enzymeIndexs_elongationFactor) && nFreTfs(this.enzymeIndexs_terminationFactor)
                trmPols{i} = actPols{i}(rnaPols.states{i}(actPols{i}) > ...
                    tuLens{i}(trnspts.boundTranscriptionUnits{i}(actPols{i})));
                nTrmPols{i} = numel(trmPols{i});
                if nTrmPols{i} > this.substrates(this.substrateIndexs_water)
                   order{i} = this.randStream.randperm(nTrmPols{i});
                   nTrmPols{i} = this.substrates(this.substrateIndexs_water);
                   trmPols{i} = trmPols{i}(order{i}(1:nTrmPols{i}));
                end
                
                if nTrmPols{i} > 0
                    if isscalar(trmPols{i})
                        this.RNAs(trnspts.boundTranscriptionUnits{i}(trmPols{i})) = ...
                            this.RNAs(trnspts.boundTranscriptionUnits{i}(trmPols{i})) + 1;
                    else
                        this.RNAs = this.RNAs + ...
                            histc(trnspts.boundTranscriptionUnits{i}(trmPols{i}, 1), 1:numel(tuLens{i}));
                    end
                    
                    %release proteins from chromsome
                    releasedProteins{i} = this.releaseProteinFromSites(i, rnaPols.positionStrands{i}(trmPols{i}, :), false, this.enzymeIndexs_rnaPolymerase1, true, true);
                    if ~isequal(releasedProteins{i}(this.enzymeIndexs_rnaPolymerase1), numel(trmPols{i}))
                        %throw(MException('Transcription:error', 'Unable to unbind protein'));
						warning('Unable to unbind protein');
                    end
					
                    %nFrePols = nFrePols + nTrmPols{i};
                    %nBndPols = nBndPols - nTrmPols{i};
                    nFrePols = nFrePols + nTrmPols{i};
                    nBndPols = nBndPols - nTrmPols{i};
                    
                    this.substrates(this.substrateIndexs_water)    = this.substrates(this.substrateIndexs_water)    - nTrmPols{i};
                    this.substrates(this.substrateIndexs_hydrogen) = this.substrates(this.substrateIndexs_hydrogen) + nTrmPols{i};
                    
                    rnaPols.states{i}(trmPols{i}) = rnaPols.freeValue;
                    rnaPols.positionStrands{i}(trmPols{i}, :) = 0;
                    trnspts.boundTranscriptionUnits{i}(trmPols{i}) = trnspts.nullTranscriptValue;
                    trnspts.boundTranscriptProgress{i}(trmPols{i}) = trnspts.nullTranscriptValue;
                    trnspts.boundTranscriptChromosome{i}(trmPols{i}) = trnspts.nullTranscriptValue;
                end
            end
            
            %% diphosphate released by polymerization
            this.substrates(this.substrateIndexs_diphosphate) = ...
                this.substrates(this.substrateIndexs_diphosphate) + ...
                sum(usedNTPs);
            end
			
            %% store enzyme state
            this.enzymes(this.enzymeIndexs_transcriptionFactors) = nFreTfs;
            this.boundEnzymes(this.enzymeIndexs_transcriptionFactors) = nBndTFs;
            this.enzymes(this.enzymeIndexs_rnaPolymerase1) = nFrePols;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase1) = nBndPols;
            this.enzymes(this.enzymeIndexs_rnaPolymerase2) = nFreHols;
            this.boundEnzymes(this.enzymeIndexs_rnaPolymerase2) = nBndHols;
		end
    end
    
    %model helper functions
    methods    
        function pBinds = computeRNAPolymeraseTUBindingProbabilities(this, i, protectDnaABoxes)
            c = this.chromosome;
            r = this.rnaPolymerases;
            
            %relative probability RNA polymerase binds each transcriptionunit
			
            pBinds = this.transcriptionUnitBindingProbabilities{i}(:, ones(1, 2)) .* ...
                r.transcriptionFactorBindingProbFoldChange{i} .* ...
                r.supercoilingBindingProbFoldChange{i};
				
            %if functional DnaA boxes R1-5 are occupied by DnaA-ATP,
            %don't allow RNA polymerase to bind since it would knockoff
            %DnaA-ATP
			
            if nargin == 2 || protectDnaABoxes
                [posStrnds, complexs] = find(...
                    c.complexBoundSites{i}(c.transcriptionUnitStartCoordinates{i}(this.transcriptionUnitIndexs_DnaAR12345Boxes{i}) + ...
                    (1:c.transcriptionUnitLengths{i}(this.transcriptionUnitIndexs_DnaAR12345Boxes{i}))' - 1, :));

                if any(ismembc(complexs, this.complexIndexs_DnaA_ATP) & posStrnds(:, 2) <= 2)
                    pBinds(this.transcriptionUnitIndexs_DnaAR12345Boxes{i}, 1) = 0;
                end
                if any(ismembc(complexs, this.complexIndexs_DnaA_ATP) & posStrnds(:, 2) > 2)
                    pBinds(this.transcriptionUnitIndexs_DnaAR12345Boxes{i}, 2) = 0;
                end									
            end
        end
        
        function pStProbs = calcStateTransitionProbabilities(this, nPols, ntpProdRate, tuBindingProbs, i, method)
            import edu.jiangnan.fmme.util.ComputationUtil;
            import edu.jiangnan.fmme.util.ConstantUtil;
			
            rna = this.rna;
            rnaPols = this.rnaPolymerases;
			c = this.chromosome;
			%for i = 1:16
			A{i}  = rna.nascentIndexs([c.genome(i).transcriptionUnits.idx]);
            pStProbs = zeros(4, 4);
			pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) = ...
                ntpProdRate / (nPols * rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex)) / ...
                (tuBindingProbs' * rna.lengths(A{i}));
            if nargin < 6 || isequal(method, 'handTuned')
                pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex) = 1;
                pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.specificallyBoundIndex) = 0;
                pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = 0.10;
            else
                pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex) = ...
                    pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) * ...
                    rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex) / ...
                    rnaPols.stateExpectations(rnaPols.specificallyBoundIndex);
                pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = ...
                    pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex) * ...
                    rnaPols.stateExpectations(rnaPols.activelyTranscribingIndex) / ...
                    rnaPols.stateExpectations(rnaPols.nonSpecificallyBoundIndex);
            end
            
            pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.activelyTranscribingIndex) = ...
                1 - pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.activelyTranscribingIndex);
            pStProbs(rnaPols.specificallyBoundIndex, rnaPols.specificallyBoundIndex) = ...
                + 1 ...
                - pStProbs(rnaPols.activelyTranscribingIndex, rnaPols.specificallyBoundIndex)...
                - pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.specificallyBoundIndex);
            pStProbs(rnaPols.nonSpecificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex) = ...
                1 - pStProbs(rnaPols.specificallyBoundIndex, rnaPols.nonSpecificallyBoundIndex);
            pStProbs(:, rnaPols.freeIndex) = pStProbs(:, rnaPols.nonSpecificallyBoundIndex);
			end
		%end
    end
    
    %get methods of dependent local state
    methods
        function value = getDryWeight(this)
            if size(this.RNAs, 3) == 1
                value = this.getDryWeight@edu.jiangnan.fmme.cell.sim.Process() + ...
                    this.rna.molecularWeights(this.rna.nascentIndexs)' * this.RNAs ...
                    / edu.jiangnan.fmme.util.ConstantUtil.nAvogadro;
            else
                value = this.getDryWeight@edu.jiangnan.fmme.cell.sim.Process() + ...
                    permute(this.rna.molecularWeights(this.rna.nascentIndexs)' * permute(this.RNAs,[1 3 2]),[1 3 2]) ...
                    / edu.jiangnan.fmme.util.ConstantUtil.nAvogadro;
            end
        end
    end
end
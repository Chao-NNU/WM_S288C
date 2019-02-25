%RNA Polymerase
%
% @wholeCellModelID State_RNAPolymerase
% @name             RNA Polymerases
% @description
%
%   states represents the current state / pseudostate (actively
%   transcribing, specifically bound, non-specifically bound, free,
%   non-existent) of each RNA polymerase, where each state is indicated by the
%   enumeration:
%   - rnaPolymeraseActivelyTranscribingValue
%   - rnaPolymeraseSpecificallyBoundValue
%   - rnaPolymeraseNonSpecificallyBoundValue
%   - rnaPolymeraseFreeValue
%   - rnaPolymeraseNotExistValue (state exists as a way to account for memory
%     allocated for future RNA polymerases)
%
% Information about positions of the polymerases on the DNA and the
% progress of RNA polymerases transcribing specific transcrips is all
% contained within the chromosomeState class and newTranscriptState class.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Affilitation: FMME Lab, School of Biotechnology, Jiangna University
% Last updated: 12/13/2010
% Last modified: 10/9/2018

classdef RNAPolymerase < edu.jiangnan.fmme.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'stateExpectations'
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames              = {   %names of properties which are part of the simulation's state
            'states'
            'positionStrands'
            'transcriptionFactorBindingProbFoldChange'
            'supercoilingBindingProbFoldChange'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'stateOccupancies'
            'nActive'
            'nSpecificallyBound'
            'nNonSpecificallyBound'
            'nFree'
            };
    end
    
    %constants
    properties (Constant)
        activelyTranscribingIndex   = 1; %index within rnaPolymeraseStateOccupancies
        specificallyBoundIndex      = 2; %index within rnaPolymeraseStateOccupancies
        nonSpecificallyBoundIndex   = 3; %index within rnaPolymeraseStateOccupancies
        freeIndex                   = 4; %index within rnaPolymeraseStateOccupancies
        
        activelyTranscribingValue   = 1; %value within states
        specificallyBoundValue      = -3; %value within states
        nonSpecificallyBoundValue   = -1; %value within states
        freeValue                   = -2; %value within states
        notExistValue               = 0; %value within states
        
        stateValues                 = [1; -3; -1; -2]; %values of states
    end
    
    %fixed biological constants
    properties
        stateExpectations       %Expected fractional occupancies of RNA polymerase states
    end
    
    %state
    properties
        states = cell(16,1)   %RNA polymerase state (see Transcription process)
        positionStrands = cell(16,1) %the position on the chromosome(s) for each RNApolymerase
        transcriptionFactorBindingProbFoldChange = cell(16,1) 
		%fold change effect of currently active transcription factors on RNA polymerase binding probabilities [nTUs X 2]
        supercoilingBindingProbFoldChange = cell(16,1)       
		%fold change in transcription probabilities due to supercoiling sigma, [# transcUnits x 2 chromosomes]
    end
    
    properties (Dependent = true)
        stateOccupancies        %number of RNA polymerases in various states
        nActive                 %number of actively transcribing RNA polymerases
        nSpecificallyBound      %number of specifically bound RNA polymerases
        nNonSpecificallyBound   %number of non-specifically bound RNA polymerases
        nFree                   %number of free RNA polymerases
    end
    
    %dependent state
    properties (Constant)
        dryWeight = 0;          %dry weight of this class' state properties
    end
    
    %references to other parts of cell state
    properties
        chromosome
        transcripts
    end
    
    %constructor
    methods
        function this = RNAPolymerase(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.chromosome = simulation.state('Chromosome');
            this.transcripts = simulation.state('Transcript');
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
		
			%this.states%{i}                                   = zeros(0, 1, numTimePoints);
			for i = 1:16
            nTUs{i} = numel(this.chromosome.transcriptionUnitWholeCellModelIDs{i});
			this.states{i}									 = zeros(0, 2, numTimePoints);
			%this.positionStrands{i}                          = zeros(nTUs{i}, 2, numTimePoints);
            this.positionStrands{i}                          = zeros(0, 2, numTimePoints);
            this.transcriptionFactorBindingProbFoldChange{i} = zeros(nTUs{i}, 2, numTimePoints);
            this.supercoilingBindingProbFoldChange{i}        = zeros(nTUs{i}, 2, numTimePoints);
			end
        end
    end
    
    %initialize state
    methods
        function initialize(~)
        end
    end
    
    %external interface
    methods
        function degradeFreePolymerase(this, nPolymerases)
			for i = 1:16
				idxs{i} = find(this.states{i} == this.freeValue, nPolymerases, 'first');
				this.states{i}(idxs{i}) = this.notExistValue;
			end
        end
        
        function releasePolymerase(this, posStrnds, proteinIsDegraded)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
			for i = 1:16
			
            c = this.chromosome;
            t = this.transcripts;
            
            tfs{i} = CircularSparseMat.ismember_subs(this.positionStrands{i}, posStrnds, [c.sequenceLen(i) c.nCompartments]);
            idxs{i} = find(tfs{i});
            if isempty(idxs{i})
                return;
            end
            
            tfs{i}(tfs{i}) = this.states{i}(tfs{i}) > this.activelyTranscribingValue;
            idxs2{i} = find(tfs{i});
            for j = 1:numel(idxs2{i})%change form i to j
                if t.boundTranscriptProgress(idxs2{i}(j)) <= 1
                    continue;
                end
                t.abortedTranscripts{i} = [
                    t.abortedTranscripts{i}
                    t.boundTranscriptionUnits{i}(idxs2{i}(j)) t.boundTranscriptProgress{i}(idxs2{i}(j))-1
                    ];
            end
            if proteinIsDegraded
                this.states{i}(idxs{i}) = this.notExistValue;
            else
                this.states{i}(idxs{i}) = this.freeValue;
            end
            this.positionStrands{i}(idxs{i}, :) = 0;
            t.boundTranscriptionUnits{i}(idxs{i}) = t.nullTranscriptValue;
            t.boundTranscriptProgress{i}(idxs{i}) = t.nullTranscriptValue;
            t.boundTranscriptChromosome{i}(idxs{i}) = t.nullTranscriptValue;
			end
        end
    end
    
    %getters
    methods
        %number of actively transcribing RNA polymerases
        function value = get.stateOccupancies(this)
            numTimePoints = size(this.states,3);
            for i = 1:16
            value{i} = zeros(4,1,numTimePoints);
            value{i}(this.activelyTranscribingIndex,:,:) = this.nActive{i};
            value{i}(this.specificallyBoundIndex,:,:)    = this.nSpecificallyBound{i};
            value{i}(this.nonSpecificallyBoundIndex,:,:) = this.nNonSpecificallyBound{i};
            value{i}(this.freeIndex,:,:)                 = this.nFree{i};
            
            value{i} = value{i} ./ repmat(sum(sum(value{i},1),2), [4,1,1]);
			end
        end
        
        %number of actively transcribing RNA polymerases
        function value = get.nActive(this)
			for i = 1:16
				value{i} = sum(this.states{i} >= this.activelyTranscribingValue);
			end
        end
        
        %number of specifically bound RNA polymerases
        function value = get.nSpecificallyBound(this)
			for i = 1:16
				value{i} = sum(this.states{i} == this.specificallyBoundValue);
			end
        end
        
        %number of non-specifically bound RNA polymerases
        function value = get.nNonSpecificallyBound(this)
			for i = 1:16
				value{i} = sum(this.states{i} == this.nonSpecificallyBoundValue);
			end
        end
        
        %number of free RNA polymerases
        function value = get.nFree(this)
			for i = 1:16
				value{i} = sum(this.states{i} == this.freeValue);
			end	
        end
    end
end

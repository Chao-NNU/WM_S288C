%Transcripts
%
% @wholeCellModelID State_Transcript
% @name             Transcripts
% @description
%
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 12/16/2010
classdef Transcript < edu.jiangnan.fmme.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity'
            'seed'
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'genomeLengths'
            'transcriptionUnitFivePrimeCoordinates'
            'transcriptionUnitDirections'
            'transcriptionUnitSequences'
            'transcriptionUnitLengths'
            };
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such  
        stateNames              = {   %names of properties which are part of the simulation's state
            'boundTranscriptionUnits'
            'boundTranscriptProgress'
            'boundTranscriptChromosome'
            'abortedTranscripts'
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'rnaBoundRNAPolymerases'
            };
    end
       
    %constants
    properties (Constant)
        nullTranscriptValue = 0; %value for transcripts that do not exist; needed to sync with rnaPolymerases
        newTranscriptValue  = 1; %value for transcripts that do not exist; needed to sync with rnaPolymerases
    end

    %constants
    properties
        genomeLengths                                 %length of genome
        
        transcriptionUnitWholeCellModelIDs            %whole cell model ids of transcription units
        transcriptionUnitFivePrimeCoordinates         %transcription unit 5' coordinates
        transcriptionUnitDirections                   %transcription unit directions
        transcriptionUnitSequences                    %transcription unit sequences
        transcriptionUnitLengths                      %transcription unit sequence lengths
        
        nmpMolecularWeights                           %molecular weights of NMPs
    end
    
    %state
    properties
        %these first three properties describing transcripts currently
        %being made need to be kept in sync with each other. 
        boundTranscriptionUnits     %transcription unit of all transcripts currently being made
        boundTranscriptProgress     %progress on transcription unit of all transcripts currently being made
        boundTranscriptChromosome   %chromosome of all transcripts currently being made
        
        abortedTranscripts          %identity of aborted transcripts (transcripts X [transcription Unit, length])
        totalBaseCounts             %total base counts of aborted transcripts
    end   
    
    %dependent state
    properties (Dependent = true)
        abortedSequences            %sequences of aborted transcripts
        rnaMaxRNAPolymeraseState    %Maximum nascent transcript length of each transcription unit
        rnaBoundRNAPolymerases      %Number of RNA polymerases actively transcribing each gene
    end

    properties (Dependent = true)
        dryWeight                   %dry weight of this class' state properties
    end
    
    %references to other parts of cell state
    properties   
        chromosome        
    end

    %constructor
    methods
        function this = Transcript(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.CellState(wholeCellModelID, name);            
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.chromosome = simulation.state('Chromosome');
        end
        
        function initializeConstants(this, knowledgeBase, simulation)
		    this.initializeConstants@edu.jiangnan.fmme.cell.sim.CellState(knowledgeBase, simulation);
			
            this.genomeLengths                          = [knowledgeBase.genomes.sequenceLength]';
			g = knowledgeBase.genomes;
			for i = 1:16
			this.transcriptionUnitWholeCellModelIDs{i}     = {g(i).transcriptionUnits.wholeCellModelID}';
			this.transcriptionUnitFivePrimeCoordinates{i}  = [g(i).transcriptionUnits.fivePrimeCoordinate]';
			this.transcriptionUnitDirections{i}            = [g(i).transcriptionUnits.direction]';
			this.transcriptionUnitSequences{i}             = {g(i).transcriptionUnits.sequence}';
			this.transcriptionUnitLengths{i}               = [g(i).transcriptionUnits.sequenceLength]';
			end			
            this.nmpMolecularWeights = simulation.state('Metabolite').molecularWeights(simulation.state('Metabolite').nmpIndexs);
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)
			for i = 1:16
            this.boundTranscriptionUnits{i}   = zeros(0, 1, numTimePoints);
            this.boundTranscriptProgress{i}   = zeros(0, 1, numTimePoints);
            this.boundTranscriptChromosome{i} = zeros(0, 1, numTimePoints);
            
            this.abortedTranscripts{i} = zeros(0, 2, numTimePoints);
			end
        end
    end
    
    methods
        function initialize(this)
			for i = 1:16
            this.abortedTranscripts{i} = zeros(0, 2);
			end
        end
    end
    
    %getters
    methods
        function value = get.dryWeight(this)
            import edu.jiangnan.fmme.util.ConstantUtil;
			for j = 1:16 
            nSeq = 0;
            seq{j} = [];
            
            %actively transcribing transcripts
            for i = 1:size(this.boundTranscriptionUnits{j}, 1)
                if this.boundTranscriptProgress{j}(i) > this.nullTranscriptValue
                    nSeq = nSeq + 1;
                    seq{j} = [seq{j} this.transcriptionUnitSequences{j}{this.boundTranscriptionUnits{j}(i)}(1:(this.boundTranscriptProgress{j}(i) - 1))]; %#ok<AGROW>
                end
            end
            
            %aborted
            for i = 1:size(this.abortedTranscripts{j}, 1)
                if this.abortedTranscripts{j}(i, 1) ~= 0
                    nSeq = nSeq + 1;
                    seq{j} = [seq{j} this.transcriptionUnitSequences{j}{this.abortedTranscripts{j}(i, 1)}(1:this.abortedTranscripts{j}(i, 2))]; %#ok<AGROW>
                end
            end
            
            %mass
            value{j} = (...
                + [sum(seq{j} == 'A') sum(seq{j} == 'C') sum(seq{j} == 'G') sum(seq{j} == 'U')] * this.nmpMolecularWeights ...
                - (numel(seq{j}) - nSeq) * (ConstantUtil.elements.H + ConstantUtil.elements.O) ...
                ) / ConstantUtil.nAvogadro;
			end
			
			value = value{1} + value{2} + value{3} + value{4} + value{5} + value{6} + value{7} + value{8} + value{9} + value{10}...
				  + value{11} + value{12} + value{13} + value{14} + value{15} + value{16};
        end
        
        function value = get.abortedSequences(this)
		
			for j = 1:16
				value{j} = cell(size(this.abortedTranscripts{j}, 1), 1);
				for i = 1:size(this.abortedTranscripts{j}, 1)
					if this.abortedTranscripts{j}(i, 1) == 0
						continue;
					end
					value{j}{i} = this.transcriptionUnitSequences{j}{this.abortedTranscripts{j}(i, 1)}(1:this.abortedTranscripts{j}(i, 2));
				end
			end
        end
        
        function value = get.totalBaseCounts(this)
            import edu.jiangnan.fmme.cell.kb.ssRNA;
            for j = 1:16
            value{j} = zeros(1, 4);
            
            for i = 1:size(this.boundTranscriptionUnits{j}, 1)
                if this.boundTranscriptProgress{j}(i) > this.nullTranscriptValue
                    seq{j} = this.transcriptionUnitSequences{j}{this.boundTranscriptionUnits{j}(i)}(1:(this.boundTranscriptProgress{j}(i) - 1));
                    value{j} = value{j} + ssRNA.computeBaseCount(seq{j}, 4, 1:4);
                end
            end
            
            abortedSeqs{j} = this.abortedSequences{j};
            for i = 1:numel(abortedSeqs{j})
                value{j} = value{j} + ssRNA.computeBaseCount(abortedSeqs{j}{i}, 4, 1:4);
            end
			end
        end
        
        %Maximum nascent transcript length of each transcription unit
        function value = get.rnaMaxRNAPolymeraseState(this)
			for j = 1:16
            value{j} = zeros(length(this.transcriptionUnitLengths{j}), 1, size(this.boundTranscriptionUnits{j}, 3));
            for i = 1:size(this.boundTranscriptionUnits{j}, 1)
                for k = 1:size(this.boundTranscriptionUnits{j}, 3)
                    boundTranscriptionUnit{j} = this.boundTranscriptionUnits{j}(i, 1, k);
                    if boundTranscriptionUnit{j} == 0
                        continue;
                    end
                    value(boundTranscriptionUnit{j}, 1, k) = max(...
                        value(boundTranscriptionUnit{j}, 1, k),...
                        this.transcriptionUnitLengths{j}(i, 1, k));
                end
            end
			end
        end

        %Number of RNA polymerases actively transcribing each gene
        function value = get.rnaBoundRNAPolymerases(this)
			for j = 1:16
            nTUs{j} = numel(this.transcriptionUnitLengths{j});
            nTimePoints{j} = size(this.boundTranscriptionUnits{j}, 3);
            value{j} = zeros(nTUs{j}, 1, nTimePoints{j});
            for k = 1:nTimePoints{j}
                value{j}(:, 1, k) = histc(this.boundTranscriptionUnits{j}(:, 1, k), 1:nTUs{j});
            end
			end
        end
    end
end

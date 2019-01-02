% Defines a parameter
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/17/2009
classdef GenomeFeature < edu.jiangnan.fmme.cell.kb.dsDNA

    properties
        genomes            = edu.jiangnan.fmme.cell.kb.Genome.empty(0, 0);
		genes              = edu.jiangnan.fmme.cell.kb.Gene.empty(0, 0);
		transcriptionUnits = edu.jiangnan.fmme.cell.kb.TranscriptionUnit.empty(0, 0);
    end

    properties %(SetAccess = protected)
        halfLife

        type
        subtype
        startCoordinate
        endCoordinate
        direction
        sequence
    end

    %computed properties
    properties %(SetAccess = protected)
        fivePrimeCoordinate
        threePrimeCoordinate
        %genes
        %transcriptionUnits
    end

    methods
        function this = GenomeFeature(knowledgeBase, wid, wholeCellModelID, name,...
                type, subtype, startCoordinate, endCoordinate, length, direction, sequence, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.jiangnan.fmme.cell.kb.GenomeFeature.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.jiangnan.fmme.cell.kb.GenomeFeature;
            for i = 1:size(wid, 1)
                this(i,1).idx = i;
                this(i,1).knowledgeBase = knowledgeBase;
                this(i,1).wid = wid(i);
                this(i,1).wholeCellModelID = wholeCellModelID{i};
                this(i,1).name = name{i};
                this(i,1).type = type{i};
                this(i,1).subtype = subtype{i};
                this(i,1).startCoordinate = startCoordinate(i);
				this(i,1).endCoordinate = endCoordinate(i);
                %this(i,1).endCoordinate = startCoordinate(i); %+ length(i) - 1;
                this(i,1).direction = direction(i);
				this(i,1).sequence = sequence{i};
                if exist('comments', 'var') && ~isempty(comments); this(i, 1).comments = comments{i}; end;
                if exist('crossReferences', 'var')
                    if size(crossReferences, 1) > 1
                        this(i, 1).crossReferences = crossReferences(i);
                    else
                        this(i, 1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields, 1)
                            values = crossReferences.(fields{j});
                            this(i, 1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genomes              = this.serializeLinksHelper(this(i).genomes);
				this(i).genes                = this.serializeLinksHelper(this(i).genes);
				this(i).transcriptionUnits   = this.serializeLinksHelper(this(i).transcriptionUnits);	
                serializeLinks@edu.jiangnan.fmme.cell.kb.dsDNA(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genomes              = this.deserializeLinksHelper(this(i).genomes, kb.genomes);
                this(i).genes                = this.deserializeLinksHelper(this(i).genes, kb.genes);
				this(i).transcriptionUnits   = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                deserializeLinks@edu.jiangnan.fmme.cell.kb.dsDNA(this(i), kb);
            end
        end

        function value = get.fivePrimeCoordinate(this)
            %retrieve
            if ~isempty(this.fivePrimeCoordinate)
                value = this.fivePrimeCoordinate;
                return;
            end
            
            %compute
            if this.direction
                value = this.startCoordinate;
            else
                value = this.endCoordinate;
            end
            
            %store
            this.fivePrimeCoordinate = value;
        end

        function value = get.threePrimeCoordinate(this)
            %retrieve
            if ~isempty(this.threePrimeCoordinate)
                value = this.threePrimeCoordinate;
                return;
            end
            
            %compute
            if this.direction
                value = this.endCoordinate;
            else
                value = this.startCoordinate;
            end
            
            %store
            this.threePrimeCoordinate = value;
        end
%{
        function value = get.sequence(this)
            %retrieve
            if ~isempty(this.sequence)
                value = this.sequence;
                return;
            end
            
            %compute
            if this.direction
                value = this.genome.sequence(this.startCoordinate:this.endCoordinate);
            else
                value = seqrcomplement(this.genome.sequence(this.startCoordinate:this.endCoordinate));
            end            
            
            %store
            this.sequence = value;
        end
        function value = get.genes(this)
            %retrieve
            if ~isempty(this.genes)
                value = this.genes;
                return;
            end
            
            %compute
            geneStartCoordinates = [this.genomes.genes.startCoordinate];
            geneEndCoordinates = [this.genomes.genes.endCoordinate];
            value = this.genomes.genes(...
                (this.startCoordinate > geneStartCoordinates & this.startCoordinate < geneEndCoordinates) | ...
                (this.endCoordinate   > geneStartCoordinates & this.endCoordinate   < geneEndCoordinates) | ...
                (this.startCoordinate < geneStartCoordinates & this.endCoordinate   > geneEndCoordinates));
                        
            %store
            this.genes = value;
        end

        function value = get.transcriptionUnits(this)
            %retrieve
            if ~isempty(this.transcriptionUnits)
                value = this.transcriptionUnits;
                return;
            end
            
            %compute
            transcriptionUnitStartCoordinates = [this.genomes.transcriptionUnits.startCoordinate];
            transcriptionUnitEndCoordinates = [this.genomes.transcriptionUnits.endCoordinate];
            value = this.genomes.transcriptionUnits(...
                (this.startCoordinate > transcriptionUnitStartCoordinates & this.startCoordinate < transcriptionUnitEndCoordinates) | ...
                (this.endCoordinate   > transcriptionUnitStartCoordinates & this.endCoordinate   < transcriptionUnitEndCoordinates) | ...
                (this.startCoordinate < transcriptionUnitStartCoordinates & this.endCoordinate   > transcriptionUnitEndCoordinates));            
            
            %store
            this.transcriptionUnits = value;
        end
%}
        %halfLife
        function value = get.halfLife(this)
            throw(MException('GenomeFeature:error', 'property is not defined'));
        end
    end
end
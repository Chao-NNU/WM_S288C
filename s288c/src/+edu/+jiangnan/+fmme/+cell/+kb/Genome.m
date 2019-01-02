% Defines a genome
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef Genome < edu.jiangnan.fmme.cell.kb.dsDNA

    properties
        genes                    = edu.jiangnan.fmme.cell.kb.Gene.empty(0, 0);
        transcriptionUnits       = edu.jiangnan.fmme.cell.kb.TranscriptionUnit.empty(0, 0);
        genomeFeatures           = edu.jiangnan.fmme.cell.kb.GenomeFeature.empty(0, 0);
    end

    properties %(SetAccess = protected)
        halfLife
        sequence
    end

    methods
        %constructor
        function this = Genome(knowledgeBase, wid, wholeCellModelID, name, ...
                sequenceTopology, sequence, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.jiangnan.fmme.cell.kb.Genome.empty(size(wid, 1), 0);
            this(size(wid, 1), 1) = edu.jiangnan.fmme.cell.kb.Genome;
            for i = 1:size(wid, 1)
                this(i, 1).idx = i;
                this(i, 1).knowledgeBase = knowledgeBase;
                this(i, 1).wid = wid(i);
                this(i, 1).wholeCellModelID = wholeCellModelID{i};
                this(i, 1).name = name{i}; 
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

                this(i, 1).sequenceTopology = sequenceTopology{i};
                this(i, 1).sequence = upper(sequence{i});
            end
        end

        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).genes              = this.serializeLinksHelper(this(i).genes);
                this(i).transcriptionUnits = this.serializeLinksHelper(this(i).transcriptionUnits);
                this(i).genomeFeatures     = this.serializeLinksHelper(this(i).genomeFeatures);

                serializeLinks@edu.jiangnan.fmme.cell.kb.dsDNA(this(i));
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                this(i).genes              = this.deserializeLinksHelper(this(i).genes, kb.genes);
                this(i).transcriptionUnits = this.deserializeLinksHelper(this(i).transcriptionUnits, kb.transcriptionUnits);
                this(i).genomeFeatures     = this.deserializeLinksHelper(this(i).genomeFeatures, kb.genomeFeatures);
                
                deserializeLinks@edu.jiangnan.fmme.cell.kb.dsDNA(this(i), kb);
            end
        end

        function value = get.halfLife(this)
            throw(MException('Genome:error', 'property is not defined'));
        end
    end
end

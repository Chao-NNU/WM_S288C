% Defines a rRNA polymer
%
% Author: Jonathan Karr
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 5/7/2009
classdef rRNA < edu.jiangnan.fmme.cell.kb.ssRNA & edu.jiangnan.fmme.cell.kb.Enzyme
    methods
        function this = rRNA(knowledgeBase, wid, wholeCellModelID, name, ...
                sequence, ...
                molecularInteraction, chemicalRegulation, subsystem, ...
                generalClassification, proteaseClassification, ...
                transporterClassification, ...
                comments, crossReferences)

            if nargin == 0; return; end;

            this = edu.jiangnan.fmme.cell.kb.rRNA.empty(size(wid,1),0);
            this(size(wid,1),1) = edu.jiangnan.fmme.cell.kb.ssRNA;
            for i = 1:size(wid,1)
                this(i,1).idx = i;
                this(i,1).knowledgeBase = knowledgeBase;
                this(i,1).wid = wid(i);
                this(i,1).wholeCellModelID = wholeCellModelID{i};
                this(i,1).name = name{i};
                if exist('comments','var') && ~isempty(comments); this(i,1).comments = comments{i}; end;
                if exist('crossReferences','var')
                    if size(crossReferences,1)>1
                        this(i,1).crossReferences = crossReferences(i);
                    else
                        this(i,1).crossReferences = struct;
                        fields = fieldnames(crossReferences);
                        for j = 1:size(fields,1)
                            values = crossReferences.(fields{j});
                            this(i,1).crossReferences.(fields{j}) = values(i);
                        end
                    end
                end

                this(i,1).sequence = sequence(i);
                this(i,1).molecularInteraction = molecularInteraction{i};
                this(i,1).chemicalRegulation = chemicalRegulation{i};
                this(i,1).subsystem = subsystem{i};
                this(i,1).generalClassification = generalClassification{i};
                this(i,1).proteaseClassification = proteaseClassification;
                this(i,1).transporterClassification = transporterClassification{i};
            end
        end
        
        function serializeLinks(this)
            for i = 1:numel(this)
                this(i).serializeLinks@edu.jiangnan.fmme.cell.kb.ssRNA();
                this(i).serializeLinks@edu.jiangnan.fmme.cell.kb.Enzyme();
            end
        end
        
        function deserializeLinks(this, kb)
            for i = 1:numel(this)
                deserializeLinks@edu.jiangnan.fmme.cell.kb.ssRNA(this(i), kb);
                deserializeLinks@edu.jiangnan.fmme.cell.kb.Enzyme(this(i), kb);
            end
        end
    end
end
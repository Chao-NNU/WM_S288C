%Transcriptional Regulation
%
% @wholeCellModelID Process_TranscriptionalRegulation
% @name             Transcriptional Regulation
% @description
%   Biology
%   ==============
%   The the rate of transcription of each transcription unit is known to be
%   regulated by proteins referred to as transcription factors which modulate
%   the affinity of the RNA polymerase for each promoter. Transcription factors
%   can both positively and negatively modulate RNA polymerase - promoter
%   affinity. Transcription factors can stabilize RNA polymerase - promoter
%   complexes by contributing additional negative free energy to the complex by
%   providing an additional surface for the RNA polymerase binding.
%   Transcription factors can destabilize RNA polymerase - promoter complexes
%   for example by sterically blocking promoters and preventing the RNA
%   polymerase from binding promoters.
%
%   This process simulates the binding of transcription factors to the promoters
%   of each transcription unit (affinity), as well as the affect of the
%   transcription factors on the recruiment of the RNA polymerase to each
%   promoter (activity). The process was built using experimentally observed fold
%   change expression effect of each transcription factor on each promoter.
%
%   Affinity
%   ++++++++++++++
%   Binds enzymes to promoters assuming:
%   1) Transcription factors have high affinity for promoters
%   2) Transcription factors bind promoters rapidly
%   3) Transcription factors bind promoters stably over time
%   4) Only 1 copy of a transcription factor can bind each promoter
%   5) Transcription factors only compete within their species for
%      promoters. That is for each transcription unit we assume each
%      transcription factor species binds a distinct promoters region.
%
%   Consequently, at each time step we simulate that each free
%   transcription factor binds randomly binds unoccupied promoters (no copy
%   of that transcription factor is already bound to the promoter). Random
%   transcription factor-promoter binding is weighted by the affinity of
%   each transcription factor for each promoter.
%
%   Because transcription factor-promoter affinities are generally not
%   experimentally observed, we base them on the transcription factor fold
%   change activities.
%
%   Activity
%   ++++++++++++++
%   The effect of bound transcription factors on the recruitment of RNA
%   polymerase and the expression of transcription units is simulated here
%   and incorporated into the calculation of the RNA polymerase
%   transcription unit promoter binding probabilities in the transcription
%   process. Specifically, the wild-type average RNA polymerase
%   transcription unit promoter binding probabilities are multiplied by the
%   binding probability fold change effects simulated in this process.
%
%   The RNA polymerase binding probability fold change is simulated for
%   each promoter as the product of the observed expression fold change
%   effects of each bound transcription factor. When a promoter is bound by
%   a single transcription factor, the net RNA polymerase binding
%   probability fold change is the observed expression fold change of that
%   transcription factor. When a promoter is bound by multiple
%   transcription units, the net RNA polymerase binding probability fold
%   change is given by the product of the individual fold change effects of
%   the bound transcription factors.
%
%   Knowledge Base
%   ==============
%   The list of transcriptional regulatory relationships is maintained in the
%   the knowledge base. As of 8/10/2010, it contained 31 such relationships
%   between 5 transcription factors and 29 transcription units containing 37
%   genes. The knowledge base was built from a variety of literature sources
%   and databases including:
%   - PUB_0096
%   - PUB_0110
%   - PUB_0112
%   - PUB_0196
%   - PUB_0418-20
%   - PUB_0433-8
%   - PUB_0505
%
%   Initialization
%   =================
%   Because we assume that transcription factors have high affinity for DNA and
%   bind DNA stably, we initialize as many transcription factors as possible to
%   the promoter-bound state. We randomly assign transcription factors to
%   promoters using their relative affinities.
%
%   Simulation
%   ==============
%   For each kind of transcription factor, bind any free transcription factors
%   to promoters that aren't already occupied by this kind of transcription
%   factor. Choose the binding sites randomly, weighted by the transcription
%   factor's affinity to them.
%
%   References
%   ==============
%   1) Lacramioara Bintu and Nicolas E Buchler and Hernan G Garcia and
%      Ulrich Gerland and Terence Hwa, Jane Kondev and Rob Phillips (2005).
%      Transcriptional regulation by the numbers: models. Curr Opin Genet
%      Dev. 15:116-24.
%   2) Lacramioara Bintu and Nicolas E Buchler and Hernan G Garcia and
%      Ulrich Gerland and Terence Hwa, Jane Kondev and Rob Phillips (2005).
%      Transcriptional regulation by the numbers: applications. Curr Opin
%      Genet Dev. 15:125-35.
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 10/19/2010
classdef TranscriptionalRegulation < edu.jiangnan.fmme.cell.sim.Process & edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
		    'tfIndexs';
            'tuIndexs';
            'tfAffinities';
            'tfActivities';
            'tfPositionStrands';
            'otherActivities';
		    };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};
        substrateWholeCellModelIDs = {};

        enzymeWholeCellModelIDs
        %enzymeIndexs_fur         %MG_236_MONOMER  ferric uptake repressor
        %enzymeIndexs_gntR        %MG_101_MONOMER  Uncharacterized HTH-type transcriptional regulator
        %enzymeIndexs_hrcA        %MG_205_DIMER    heat-inducible transcription repressor HrcA, putative
        %enzymeIndexs_luxR        %MG_428_DIMER    LuxR bacterial regulatory protein, putative
        %enzymeIndexs_spx         %MG_127_MONOMER  Spx subfamily protein
		
		enzymeIndexs_HSF1		  	%YGL073W_MONOMER
		enzymeIndexs_ABF1			%YKL112W_MONOMER	
		enzymeIndexs_RAP1			%YNL216W_MONOMER
		enzymeIndexs_MCM1			%YMR043W_MONOMER
		enzymeIndexs_PDR1			%YGL013C_MONOMER
		enzymeIndexs_PDR3			%YBL005W_MONOMER
		enzymeIndexs_GAL4			%YPL248C_MONOMER
		enzymeIndexs_MIG1			%YGL035C_MONOMER
		enzymeIndexs_REB1			%YBR049C_MONOMER
		enzymeIndexs_BAS1			%YKR099W_MONOMER
		enzymeIndexs_PHO2			%YDL106C_MONOMER
		enzymeIndexs_PHO4			%YFR034C_MONOMER
		enzymeIndexs_SPT15			%YER148W_MONOMER	
		enzymeIndexs_LEU3			%YLR451W_MONOMER
		enzymeIndexs_STE12			%YHR084W_MONOMER
		enzymeIndexs_CBF1			%YJR060W_MONOMER
		%enzymeIndexs_SWI6			%YLR182W_MONOMER
		enzymeIndexs_SWI5			%YDR146C_MONOMER
		enzymeIndexs_FKH2			%YNL068C_MONOMER
		enzymeIndexs_ADR1			%YDR216W_MONOMER
		enzymeIndexs_CAT8			%YMR280C_MONOMER
		enzymeIndexs_HAC1			%YFL031W_MONOMER
		enzymeIndexs_NDT80			%YHR124W_MONOMER
		enzymeIndexs_UME6			%YDR207C_MONOMER
		enzymeIndexs_PUT3			%YKL015W_MONOMER
		enzymeIndexs_CUP2			%YGL166W_MONOMER
		enzymeIndexs_SKO1			%YNL167C_MONOMER
		enzymeIndexs_DAL82			%YNL314W_MONOMER
		enzymeIndexs_RME1			%YGR044C_MONOMER
		enzymeIndexs_ACE2			%YLR131C_MONOMER
		enzymeIndexs_RFA1			%YAR007C_YNL312W_YJL173C_TRIMER
		enzymeIndexs_HAP			%YGL237C_YBL021C_YKL109W_YOR358W_TETRAMER
		enzymeIndexs_INO			%YDR123C_YOL108C_DIMER
		enzymeIndexs_SWI			%YER111C_YLR182W_DIMER


        transcriptionUnitWholeCellModelIDs
    end

    %fixed biological constants
    properties
        tfIndexs           %identifying index of transcription factor that binds each site [sites X 2]
        tuIndexs           %identifying index of transcription unit associated with each TF binding site [sites X 2]
        tfAffinities       %affinity of each transcription factor for each promoter it binds [sites X 2]
        tfActivities       %fold change expression of each transcription unit when its promoter is bound by a particular transcription factor [sites X 2]
        tfPositionStrands  %start position and strand of each transcription factor DNA binding site [2*sites X 2]
        otherActivities    %fold change each non-DNA-binding transcription factor has on each transcription unit's expression [nTFs X nTUs]
    end

    %dependent local state (implemented as dependent property for convenience)
    properties (Dependent, SetAccess = protected)
        tfBoundPromoters       %number of each transcription factor bound to promoter of each transcription unit [transcription factors X promoters X 2]
        boundTFs               %count of how many of each transcription factor are bound [transcription factors X 1]
    end
    
    %references to cell state
    properties
        rnaPolymerase             %reference to RNA polymerase state
    end

    %constructor
    methods
        function this = TranscriptionalRegulation(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.Process(wholeCellModelID, name);
        end
    end

    %communication between process/simulation
    methods
        %set references to state objects
        function storeObjectReferences(this, simulation)
            this.storeObjectReferences@edu.jiangnan.fmme.cell.sim.Process(simulation);
            this.storeObjectReferences@edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect(simulation);
            
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.states = [this.states; {this.rnaPolymerase}];
        end
        
        %initialize constants
        function initializeConstants(this, knowledgeBase, varargin)
            %build object arrays of transcription factor monomers and complexes
			monomers = cell(17,1);
			complexs = cell(17,1);
			for i = 1:length(knowledgeBase.genomes)
				m{i} = {knowledgeBase.genomes(i).transcriptionUnits.transcriptionFactorProteinMonomers}';
				c{i} = {knowledgeBase.genomes(i).transcriptionUnits.transcriptionFactorProteinComplexs}';
				for j = 1:length(m{i})
					if ~isempty(m{i}{j})
						monomers{i} = [monomers{i}; m{i}{j}]; %#ok<AGROW>
					end
				end
				for k = 1:length(c{i})
					if ~isempty(c{i}{k})
						complexs{i} = [complexs{i}; c{i}{k}]; %#ok<AGROW>
				end
				end
			end
			
			monomer = [];
            complex = [];
					
            for i = 1:knowledgeBase.numTranscriptionUnits
                mo = knowledgeBase.transcriptionUnits(i).transcriptionFactorProteinMonomers;
                co = knowledgeBase.transcriptionUnits(i).transcriptionFactorProteinComplexs;
                if ~isempty(mo)
                    monomer = [monomer; mo]; %#ok<AGROW>
                end
                if ~isempty(co)
                    complex = [complex; co]; %#ok<AGROW>
                end
            end
            %whole cell model ids of all transcription factors
            this.enzymeWholeCellModelIDs = unique({
                monomer.wholeCellModelID...
                complex.wholeCellModelID}');				


            %transcription factors indices
			
			this.enzymeIndexs_HSF1  = this.enzymeIndexs({'YGL073W_MONOMER'});
			this.enzymeIndexs_ABF1  = this.enzymeIndexs({'YKL112W_MONOMER'});
			this.enzymeIndexs_RAP1  = this.enzymeIndexs({'YNL216W_MONOMER'});
			this.enzymeIndexs_MCM1  = this.enzymeIndexs({'YMR043W_MONOMER'});
			this.enzymeIndexs_PDR1  = this.enzymeIndexs({'YGL013C_MONOMER'});
			this.enzymeIndexs_PDR3  = this.enzymeIndexs({'YBL005W_MONOMER'});
			this.enzymeIndexs_GAL4  = this.enzymeIndexs({'YPL248C_MONOMER'});
			this.enzymeIndexs_MIG1  = this.enzymeIndexs({'YGL035C_MONOMER'});
			this.enzymeIndexs_REB1  = this.enzymeIndexs({'YBR049C_MONOMER'});
			this.enzymeIndexs_BAS1  = this.enzymeIndexs({'YKR099W_MONOMER'});
			this.enzymeIndexs_PHO2  = this.enzymeIndexs({'YDL106C_MONOMER'});
			this.enzymeIndexs_PHO4  = this.enzymeIndexs({'YFR034C_MONOMER'});
			this.enzymeIndexs_SPT15  = this.enzymeIndexs({'YER148W_MONOMER'});
			this.enzymeIndexs_LEU3  = this.enzymeIndexs({'YLR451W_MONOMER'});
			this.enzymeIndexs_STE12  = this.enzymeIndexs({'YHR084W_MONOMER'});
			this.enzymeIndexs_CBF1  = this.enzymeIndexs({'YJR060W_MONOMER'});
			%this.enzymeIndexs_SWI6  = this.enzymeIndexs({'YLR182W_MONOMER'});
			this.enzymeIndexs_SWI5  = this.enzymeIndexs({'YDR146C_MONOMER'});
			this.enzymeIndexs_FKH2  = this.enzymeIndexs({'YNL068C_MONOMER'});
			this.enzymeIndexs_ADR1  = this.enzymeIndexs({'YDR216W_MONOMER'});
			this.enzymeIndexs_CAT8  = this.enzymeIndexs({'YMR280C_MONOMER'});
			this.enzymeIndexs_HAC1  = this.enzymeIndexs({'YFL031W_MONOMER'});
			this.enzymeIndexs_NDT80  = this.enzymeIndexs({'YHR124W_MONOMER'});
			this.enzymeIndexs_UME6  = this.enzymeIndexs({'YDR207C_MONOMER'});
			this.enzymeIndexs_PUT3  = this.enzymeIndexs({'YKL015W_MONOMER'});
			this.enzymeIndexs_CUP2  = this.enzymeIndexs({'YGL166W_MONOMER'});
			this.enzymeIndexs_SKO1  = this.enzymeIndexs({'YNL167C_MONOMER'});
			this.enzymeIndexs_DAL82  = this.enzymeIndexs({'YNL314W_MONOMER'});
			this.enzymeIndexs_RME1  = this.enzymeIndexs({'YGR044C_MONOMER'});
			this.enzymeIndexs_ACE2  = this.enzymeIndexs({'YLR131C_MONOMER'});
			this.enzymeIndexs_RFA1  = this.enzymeIndexs({'YAR007C_YNL312W_YJL173C_TRIMER'});
			this.enzymeIndexs_HAP  = this.enzymeIndexs({'YGL237C_YBL021C_YKL109W_YOR358W_TETRAMER'});
			this.enzymeIndexs_INO  = this.enzymeIndexs({'YDR123C_YOL108C_DIMER'});
			this.enzymeIndexs_SWI  = this.enzymeIndexs({'YER111C_YLR182W_DIMER'});

            %call super class method to compute mapping between enzymes of
            %this process, and of protein monomers and complexes of the simulation
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.Process(...
                knowledgeBase, varargin{:});
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, varargin{:});
			c = this.chromosome;
			
            %whole cell model ids of transcription units
			
			this.tuIndexs			 = cell(16,1);
			this.tfIndexs 			 = cell(16,1);
			this.tfActivities  		 = cell(16,1);
			this.tfPositionStrands 	 = cell(16,1);
			this.otherActivities	 = cell(16,1);
           for i = 1:length(knowledgeBase.genomes)
				if ~isempty(monomers{i})|~isempty(complexs{i})
				n(i) = length(monomers{i}) + length(complexs{i});
				this.tuIndexs{i} = zeros(n(i), 2);
				this.tfIndexs{i} = zeros(n(i), 2);
				this.tfActivities{i} = zeros(n(i), 2);
				this.tfPositionStrands{i} = zeros(n(i)*2, 2);
				
				tu = knowledgeBase.genomes(i).transcriptionUnits;
				g  = tu.genomes;
				this.transcriptionUnitWholeCellModelIDs{i} = {tu.wholeCellModelID};
				j = 1;
				for z = 1:length(tu)
					if ~isempty(tu(z).transcriptionFactorProteinMonomers)
					[tf{i}, proteinIndexs{i}] = ismember(...
						[tu(z).transcriptionFactorProteinMonomers.idx]',...
						this.enzymeMonomerGlobalIndexs);
					k = nnz(tf{i}) - 1;
					this.tuIndexs{i}(j:j+k) = z;
					proteinIndexs{i} = proteinIndexs{tf{i}};
					this.tfIndexs{i}(j:j+k) = this.enzymeMonomerLocalIndexs(proteinIndexs{i});
					
					regulationIndexs{i} = ...
						this.enzymeMonomerCompartmentIndexs(proteinIndexs{i}) == ...
						[tu(z).transcriptionFactorProteinMonomerCompartments.idx]';
					this.tfActivities{i}(j:j+k) = ...
						tu(z).transcriptionFactorProteinMonomerActivitys(regulationIndexs{i});
					
					footprints{i} = c.monomerDNAFootprints([tu(z).transcriptionFactorProteinMonomers.idx]');	
					startCoordinates{i} = tu(z).transcriptionFactorProteinMonomerBindingSiteStartCoordinates;
					lengths{i} = tu(z).transcriptionFactorProteinMonomerBindingSiteLengths;
					this.tfPositionStrands{i}(j:j+k) = ...
						mod(ceil(startCoordinates{i} + lengths{i}/2 - footprints{i}/2) - 1, g.sequenceLength) + 1;
					j = j + 1 + k;
					end
					
					
					if ~isempty(tu(z).transcriptionFactorProteinComplexs)
					[tf{i}, proteinIndexs{i}] = ismember(...
						[tu(z).transcriptionFactorProteinComplexs.idx]',...
						 this.enzymeComplexGlobalIndexs);
					
					k = nnz(tf{i}) - 1;
					this.tuIndexs{i}(j:j+k) = z;
					%proteinIndexs{i} = proteinIndexs{tf{i}};
					
					this.tfIndexs{i}(j:j+k) = this.enzymeComplexLocalIndexs(proteinIndexs{i});
					
					regulationIndexs{i} = ...
						this.enzymeComplexCompartmentIndexs(proteinIndexs{i}) == ...
						[tu(z).transcriptionFactorProteinComplexCompartments.idx]';
						
					this.tfActivities{i}(j:j+k) = ...
						tu(z).transcriptionFactorProteinComplexActivitys(regulationIndexs{i});
					
					footprints{i} = c.complexDNAFootprints([tu(z).transcriptionFactorProteinComplexs.idx]');	
					startCoordinates{i} = tu(z).transcriptionFactorProteinComplexBindingSiteStartCoordinates;
					lengths{i} = tu(z).transcriptionFactorProteinComplexBindingSiteLengths;
					this.tfPositionStrands{i}(j:j+k) = ...	
						mod(ceil(startCoordinates{i} + lengths{i}/2 - footprints{i}/2) - 1, g.sequenceLength) + 1;
					j = j + 1 + k;	
					end
				end
				
				validateattributes(this.tfActivities{i}, {'numeric'}, {'nonnegative'});
				
				this.tuIndexs{i}(:,2) = this.tuIndexs{i}(:,1);
				this.tfIndexs{i}(:,2) = this.tfIndexs{i}(:,1);
				this.tfActivities{i}(:,2) = this.tfActivities{i}(:,1);

				this.tfPositionStrands{i}(1:n(i),2) = 1;
				this.tfPositionStrands{i}(n(i)+1:2*n(i),1) = this.tfPositionStrands{i}(1:n(i),1);
				this.tfPositionStrands{i}(n(i)+1:2*n(i),2) = 3;	
				
				this.tfAffinities{i} = this.tfActivities{i};
				this.tfActivities{i}(this.tfActivities{i}<1 & this.tfAffinities{i}>0) = ...
                this.tfActivities{i}(this.tfActivities{i}<1 & this.tfAffinities{i}>0) .^ -1;
				
				bind{i} = ~isnan(this.tfPositionStrands{i}(1:n(i),1));
				this.otherActivities{i} = ones(...
					length(this.transcriptionUnitWholeCellModelIDs{i}), ...
					length(this.enzymeWholeCellModelIDs));
					
				this.otherActivities{i}(sub2ind(...
					size(this.otherActivities{i}), ...
					this.tuIndexs{i}(~bind{i}), ...
					this.tfIndexs{i}(~bind{i}))) = ...
					this.tfActivities{i}(~bind{i});
				this.tuIndexs{i} = this.tuIndexs{i}(bind{i},:);
				this.tfIndexs{i} = this.tfIndexs{i}(bind{i},:);
				this.tfAffinities{i} = this.tfAffinities{i}(bind{i},:);
				this.tfActivities{i} = this.tfActivities{i}(bind{i},:);
				this.tfPositionStrands{i} = this.tfPositionStrands{i}(repmat(bind{i}, [2 1]), :);
				
			end	
			end
		end
	end

    %model
    methods
        %Calculate
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, ~, ~)
            %no substrate and byproducts required
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            
            %no enzymes required
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
        end

        %initialization:
        %- Simulation initializeState method initializes 1 undamaged chromosome
        %- Various processes may bind proteins to chromosome
        %- Here we bind transcription factors to the initial chromosome
        function initializeState(this)
            this.evolveState();
        end

        %resource requirements
        function result = calcResourceRequirements_Current(this)
            result = zeros(size(this.substrates));
        end

        %simulation
        function evolveState(this)
            %Bind transcription factors to accessible sites
			for i =1:16
            tfBoundPromoters{i} = this.bindTranscriptionFactors(i);
			
            %Compute new fold changes
			this.rnaPolymerase.transcriptionFactorBindingProbFoldChange{i} = ...
                this.calcBindingProbabilityFoldChange(i, tfBoundPromoters{i});
			end
        end
        
        %Bind transcription factors to accessible sites.
        function tfBoundPromoters = bindTranscriptionFactors(this, j)
            %Estimate which binding sites exist and are accessible to the
            %transcription unit. (Note chromosome class which ensure that sites
            %are truly accessible, and only allow transcription factors to bind
            %to these sites).
            tfBoundPromoters = this.tfBoundPromoters{j};
			
            accessible = ~tfBoundPromoters & reshape(this.chromosome.isRegionPolymerized(j, this.tfPositionStrands{j}, 1, false), [], 2);
			
			for i = 1:size(this.enzymes, 1)%need to think about the TFs differences
                %For each transcription factor:
                %  Stochastically pick transcription unit promoters elements to
                %  bind from among free promoter elements weighted by the
                %  transcription factors' affinities for the promoters. (Note:
                %  each transcription factor is assumed to bind a distinct
                %  promoter region near each transcription unit such that
                %  transcription factors only compete within their species
                %  for binding to promoters.)
                if this.enzymes(i) <= 0
                    continue;
                end
				
                sites = find((this.tfIndexs{j} == i) & accessible);
				
                if isempty(sites)
                    continue;
                end
				
                tf = this.bindProteinToChromosome(j,...
                    this.tfPositionStrands{j}(sites, :), i, this.enzymes(i), this.tfAffinities{j}(sites), [], false, 1, false, []);
				
                tfBoundPromoters(sites(tf)) = true;
			end
        end
		
		function foldChange = calcBindingProbabilityFoldChange(this, i, tfBoundPromoters)
		
            %reset fold change
			
            foldChange = ones(size(this.chromosome.transcriptionUnitStartCoordinates{i}, 1), 2);
            %fold change effect of bound TFs
            for k = 1:size(this.tfIndexs{i}, 1)
                for j = 1:size(this.tfIndexs{i}, 2)
                    if ~tfBoundPromoters(k, j)
                        continue;
                    end
                    foldChange(this.tuIndexs{i}(k, j), j) = foldChange(this.tuIndexs{i}(k, j), j) * this.tfActivities{i}(k, j);
                end
            end
            otherFoldChanges = prod(this.otherActivities{i} .^ (this.enzymes(:, ones(1, size(foldChange, 1)))' > 0), 2);
            foldChange = foldChange .* [otherFoldChanges otherFoldChanges];
        end
    end
    
    %get methods of dependent local state
    methods
        function result = get.tfBoundPromoters(this)
			for i = 1:16
			result{i} = reshape(this.isDnaBound(i, this.tfPositionStrands{i}, this.tfIndexs{i}(:)), [], 2);
			end
        end
		
        function result = get.boundTFs(this)
			result = cell(16,1);
			for i = 1:16
            result{i} = histc(this.tfIndexs{i}(logical(this.tfBoundPromoters{i}(:,1))), 1:numel(this.enzymes));
			end
		end
    end
end
%Chromosome
% Integration point for processes which interact with specific
% positions/strands of the cell's chromosome(s).
% - Represents the portion of chromosome(s) accessible to enzymes. That is
%   positions/strands which are NOT
%   - damaged in any way (no gap sites, abasic sites, damaged
%     sugar-phosphates, damaged bases, cross links, strand breaks,
%     or Holliday junctions)
%   - stably bound by enzymes
%   - single stranded
%
% Terminology:
% ==================
%         Site  single base/bond of chromosomes, indicated by strand index and
%               number of bases/bonds along 5'->3' strand from ORI [position X
%               strand]
%       Region  contiguous set of bases/bonds of chromosomes, indicated by start
%               and end positions (bases/bonds along 5'->3' strand from ORI and
%               strand (positive/negative)
%   Accessible  polymerized, not bound by protein, and not damaged
% Inaccessible  not polymerized, bound by protein, or damaged
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010

%TODO
%- more precise ignoreDamageFilter
%- include gapsites in isRegionPolymerized
classdef Chromosome < edu.jiangnan.fmme.cell.sim.CellState
    %constants
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
            'doubleStrandBreakSeparation';
            'strandBreakClassification_doubleStrandBreakSeparation';
            'strandBreakClassification_segmentLength';
            'sequence';
            'sequenceLen';
			'ARSs';
			'ARSStartpositions';
			'ARSTerpositions';
            %'oriCPosition';
            %'terCPosition';
            'transcriptionUnitStartCoordinates';
            'transcriptionUnitLengths';
            'transcriptionUnitStrands';
            'monomerDNAFootprints';
            'complexDNAFootprints';
            'monomerDNAFootprintBindingStrandedness';
            'complexDNAFootprintBindingStrandedness';
            'monomerDNAFootprintRegionStrandedness';
            'complexDNAFootprintRegionStrandedness';
            'reactionBoundMonomer';
            'reactionBoundComplex';
            'reactionMonomerCatalysisMatrix';
            'reactionComplexCatalysisMatrix';
            'reactionThresholds';
            'relaxedBasesPerTurn';
            'equilibriumSuperhelicalDensity';
            'supercoiledSuperhelicalDensityTolerance';
            };
        fittedConstantNames = {};  %names of process properties that are considered fitted constants, and should be stored with the simulation as such
        stateNames = {             %names of properties which are part of the simulation's state
            'polymerizedRegions';
            'linkingNumbers';
            'monomerBoundSites';
            'complexBoundSites';
            'gapSites';
            'abasicSites';
            'damagedSugarPhosphates';
            'damagedBases';
            'intrastrandCrossLinks';
            'strandBreaks';
            'hollidayJunctions';
            'segregated';
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            %'unpolymerizedRegions'
            'singleStrandedRegions'
            'doubleStrandedRegions'
            'geneCopyNumbers'
            %'polymerizedGenes'
            'transcriptionUnitCopyNumbers'
            %'polymerizedTranscriptionUnits'
            %'geneCopyNumbers_Accessible'
            %'transcriptionUnitCopyNumbers_Accessible'
            %'accessibleGenes'
            %'accessibleTranscriptionUnits'
            'ploidy'
            %'linkingNumbers_minFreeEnergy'
            %'supercoils'
            'superhelicalDensity'
            %'supercoiled'
            %'damagedSites'
            %'damagedSites_shifted_incm6AD'
            %'damagedSites_nonRedundant'
            %'damagedSites_excm6AD'
            %'gapSites3'
            %'gapSites5'
            %'abasicSites3'
            %'abasicSites5'
            %'damagedSugarPhosphates3'
            %'damagedSugarPhosphates5'
            %'damagedBases3'
            %'damagedBases5'
            %'strandBreaks3'
            %'strandBreaks5'
            %'intrastrandCrossLinks3'
            %'intrastrandCrossLinks5'
            %'hollidayJunctions3'
            %'hollidayJunctions5'
            %'singleStrandBreaks'
            %'doubleStrandBreaks'
            %'strandBreakClassification'
            %'restrictableMunIRMSites'
            %'hemiunmethylatedMunIRMSites'
            };
        
        nCompartments = 4;              %number of strands
        nChromosomes = 2;               %number of chromosomes
        
        strandIndexs_positive    = [1; 3]; %5'->3' strands
        strandIndexs_negative    = [2; 4]; %3'->5' strands
        
        strandIndexs_ch1         = [1; 2]; %chromosome 1 strands; sort to "mother" cell during division
        strandIndexs_ch2         = [3; 4]; %chromosome 2 strands; sort to "daughter" cell during division
        
        strandIndexs_template    = [1; 4]; %strands which are serving as templates for replication
        strandIndexs_nonTemplate = [2; 3]; %strands which are not serving as templates for replication
        
        strandIndexs_old         = [1; 4]; %strands which (at end the end of the cell cycle) were present at the beginning of the cell cycle
        strandIndexs_new         = [2; 3]; %strands which (at end the end of the cell cycle) were synthesized during the cell cycle
        
        dnaStrandedness_ssDNA = 1; %values of ssDNA footprints within *DNAFootprintStrandedness
        dnaStrandedness_dsDNA = 2; %values of dsDNA footprints within *DNAFootprintStrandedness
        dnaStrandedness_xsDNA = 3; %values of ssDNA/dsDNA footprints within *DNAFootprintStrandedness
        
        strandBreakClassification_index_NB    = 1; %index within strandBreakClassification of NB
        strandBreakClassification_index_SSB   = 2; %index within strandBreakClassification of SSB
        strandBreakClassification_index_SSB_  = 3; %index within strandBreakClassification of SSB+
        strandBreakClassification_index_2SSB  = 4; %index within strandBreakClassification of 2SSB
        strandBreakClassification_index_DSB   = 5; %index within strandBreakClassification of DSB
        strandBreakClassification_index_DSB_  = 6; %index within strandBreakClassification of DSB+
        strandBreakClassification_index_DSB__ = 7; %index within strandBreakClassification of DSB++
    end
    
    %computed ids, names, indices
    properties
        transcriptionUnitWholeCellModelIDs     %whole cell model IDs of transcription units
        transcriptionUnitNames                 %names of transcription units
        
        %monomerIndexs_ligase                   %index within ProteinMonomer.matureIndexs of DNA ligase
        %complexIndexs_dnaPolymerase            %index within ProteinComplex.matureIndexs of DNA polymerase
        %complexIndexs_DisA                     %index within ProteinComplex.matureIndexs of DisA
        %complexIndexs_rnaPolymerase            %index within ProteinComplex.matureIndexs of RNA polymerase
		
		complexIndexs_ligase4					%index within ProteinComplex.matureIndexs of DNA ligase 4
		monomerIndexs_ligase1					%index within ProteinMonomer.matureIndexs of DNA ligase 1
		complexIndexs_dnaPolA					%index within ProteinComplex.matureIndexs of DNA polymerase Alpha
		complexIndexs_dnaPolD					%index within ProteinComplex.matureIndexs of DNA polymerase Delta
		complexIndexs_dnaPolE					%index within ProteinComplex.matureIndexs of DNA polymerase Epsilon
		monomerIndexs_dnaPolG					%index within ProteinMonomer.matureIndexs of DNA polymerase Gamma
		complexIndexs_dnaPolZ					%index within ProteinComplex.matureIndexs of DNA polymerase Zeta
		monomerIndexs_BRE1						%index within ProteinMonomer.matureIndexs of BRE1
		monomerIndexs_RAD9						%index within ProteinMonomer.matureIndexs of RAD9
		monomerIndexs_RAD6           			%index within ProteinMonomer.matureIndexs of RAD6
		monomerIndexs_CHK1						%index within ProteinMonomer.matureIndexs of CHK1
		monomerIndexs_MUS81						%index within ProteinMonomer.matureIndexs of MUS81	
		monomerIndexs_DPB11						%index within ProteinMonomer.matureIndexs of DPB11
		monomerIndexs_CHS5						%index within ProteinMonomer.matureIndexs of CHS5
		complexIndexs_rnaPol1					%index within ProteinComplex.matureIndexs of RNA polymerase I
		complexIndexs_rnaPol2					%index within ProteinComplex.matureIndexs of RNA polymerase II
		complexIndexs_rnaPol3					%index within ProteinComplex.matureIndexs of RNA polymerase III
		
		
        monomerIndexs_replisome                %index within ProteinMonomer.matureIndexs of replication machinery
        complexIndexs_replisome                %index within ProteinComplex.matureIndexs of replication machinery
		
        
        reactionWholeCellModelIDs              %IDs of bound protein release reactions
        reactionNames                          %names of bound protein release reactions
    end
    
    %constants
    properties
        doubleStrandBreakSeparation                           = 1;   %max separtion in bases between SSB's to be consider a DBSB
        strandBreakClassification_doubleStrandBreakSeparation = 10;  %maximum separation of single strand breaks which is considered double strand break [PUB_0486]
        strandBreakClassification_segmentLength               = 216; %length of segments for which SSB/SSB+/2SSB/DSB/DSB+/DSB++ classification applies [PUB_0486]        
        
        sequence                               %chromosome sequence
        sequenceLen                            %chromosome sequence length (number of bases)
        sequenceGCContent                      %chromosome G/C content
		ARSs								   %Autonomously replicating sequence
        ARSStartpositions                      %Autonomously replicating sequence startCoordinates
		ARSTerpositions		                   %Autonomously replicating sequence endCoordinates
		%oriCPosition                          %oriC position
        %terCPosition                          %terC position
        
        transcriptionUnitStartCoordinates       %genomic coordinates of transcription units
        transcriptionUnitLengths                %genomic coordinates of transcription units
        transcriptionUnitStrands                %genomic direction of transcription units
        
        monomerDNAFootprints                    %number of bases of DNA each monomer occupies when bound to DNA
        complexDNAFootprints                    %number of bases of DNA each complex occupies when bound to DNA
        monomerDNAFootprintBindingStrandedness  %enumeration of binding strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        complexDNAFootprintBindingStrandedness  %enumeration of binding strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        monomerDNAFootprintRegionStrandedness   %enumeration of region strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        complexDNAFootprintRegionStrandedness   %enumeration of region strandedness of DNA footprints (ssDNA, dsDNA, ssDNA/dsDNA) using this.dnaStrandedness_*sDNA
        
        reactionBoundMonomer                    %protein monomer released by each reaction [reactions X 1]
        reactionBoundComplex                    %protein complex released by each reaction [reactions X 1]
        reactionMonomerCatalysisMatrix          %monomers required to catalyze the release of a bound protein [reactions X monomers]
        reactionComplexCatalysisMatrix          %complexs required to catalyze the release of a bound protein [reactions X complexs]
        reactionThresholds                      %number of proteins required to catalyze each release reaction [reactions X 1]
        
        relaxedBasesPerTurn                     %Number of dna bases per turn for the relaxed LK calculation (10.5)
        equilibriumSuperhelicalDensity          %equilibrium superhelical density; also known as specific linking difference / \sigma_{sp} [-0.06; PUB_0749]
        supercoiledSuperhelicalDensityTolerance %tolerance in superhelical density to be considered supercoiled (0.1)
    end
    
    %state
    properties
        polymerizedRegions 			= cell(16,1)  %integers [positions x strands] indicating the start positions of polymerized regions of strands and their lengths
        linkingNumbers       		= cell(16,1)  %integers [positions x strands] indicating the current linking number of each double-stranded region
        
        monomerBoundSites        	= cell(16,1)  %indices [positions x strands] indicating start positions of protein monomers bound to DNA bases
        complexBoundSites        	= cell(16,1)  %indices [positions x strands] indicating start positions of macromolecular complexes bound to DNA bases
        
        gapSites                    = cell(16,1)  %boolean [positions x strands] indicating positions of gap sites
        abasicSites              	= cell(16,1)  %boolean [positions x strands] indicating positions of abasic sites
        damagedSugarPhosphates      = cell(16,1)  %indices [positions x strands] indicating metabolite identity of damaged sugar-phosphates
        damagedBases                = cell(16,1)  %indices [positions x strands] indicating metabolite identity of damaged bases
        intrastrandCrossLinks       = cell(16,1)  %boolean [positions x strands] indicating metabolite identity of intrastrand cross links in DNA
        strandBreaks                = cell(16,1)  %boolean [positions x strands] indicating positions of strand breaks in strands of DNA
        hollidayJunctions           = cell(16,1)  %boolean [positions x strands] indicating positions of holliday junctions
        
        segregated                  = zeros(16,1) %boolean indicating whether or not the chromsomes are segregated
    end
    
    %properties to keep track of whether or not the dependent properties need to
    %be recomputed
    properties (SetAccess = protected)
		validated = ones(16,1);
        
        validated_polymerizedRegions 						= cell(16,1)
        validated_linkingNumbers	 						= cell(16,1)
        validated_proteinBoundSites							= cell(16,1)
        validated_damaged									= cell(16,1)
        validated_gapSites									= cell(16,1)
        validated_abasicSites								= cell(16,1)
        validated_damagedSugarPhosphates					= cell(16,1)
        validated_damagedBases								= cell(16,1)
        validated_intrastrandCrossLinks						= cell(16,1)
        validated_strandBreaks								= cell(16,1)
        validated_hollidayJunctions							= cell(16,1)
        validated_segregated								= cell(16,1)
        
        validated_unpolymerizedRegions						= cell(16,1)
        validated_singleStrandedRegions						= cell(16,1)
        validated_doubleStrandedRegions						= cell(16,1)
        validated_geneCopyNumbers							= cell(16,1)
        validated_ploidy									= cell(16,1)
        validated_polymerizedGenes							= cell(16,1)
        validated_transcriptionUnitCopyNumbers				= cell(16,1)
        validated_polymerizedTranscriptionUnits				= cell(16,1)
        validated_geneCopyNumbers_Accessible				= cell(16,1)
        validated_transcriptionUnitCopyNumbers_Accessible	= cell(16,1)
        validated_accessibleGenes							= cell(16,1)
        validated_accessibleTranscriptionUnits				= cell(16,1)
        validated_linkingNumbers_minFreeEnergy				= cell(16,1)
        validated_supercoils								= cell(16,1)
        validated_superhelicalDensity						= cell(16,1)
        validated_supercoiled								= cell(16,1)
        validated_damagedSites								= cell(16,1)	
        validated_damagedSites_shifted_incm6AD				= cell(16,1)
        validated_damagedSites_nonRedundant					= cell(16,1)
        validated_damagedSites_excm6AD						= cell(16,1)
        validated_gapSites3									= cell(16,1)
        validated_gapSites5									= cell(16,1)
        validated_abasicSites3								= cell(16,1)
        validated_abasicSites5								= cell(16,1)			
        validated_damagedSugarPhosphates3					= cell(16,1)				
        validated_damagedSugarPhosphates5					= cell(16,1)
        validated_damagedBases3								= cell(16,1)
        validated_damagedBases5								= cell(16,1)
        validated_strandBreaks3								= cell(16,1)
        validated_strandBreaks5								= cell(16,1)
        validated_intrastrandCrossLinks3					= cell(16,1)
        validated_intrastrandCrossLinks5					= cell(16,1)
        validated_hollidayJunctions3						= cell(16,1)
        validated_hollidayJunctions5						= cell(16,1)
        validated_singleStrandBreaks						= cell(16,1)	
        validated_doubleStrandBreaks						= cell(16,1)
        validated_strandBreakClassification					= cell(16,1)
        validated_munIRMSiteMethylationStatus				= cell(16,1)
        validated_munIRMSiteRestrictionStatus				= cell(16,1)	
        validated_dryWeight									= cell(16,1)
    end
    
    %dependent local state (implemented as dependent property for convenience)
    properties (SetAccess = protected)
        unpolymerizedRegions %integers indicating the start positions of unpolymerized regions (ie. not yet replicated) of strands and their lengths
        singleStrandedRegions 
        doubleStrandedRegions
        geneCopyNumbers
        polymerizedGenes
        transcriptionUnitCopyNumbers
        polymerizedTranscriptionUnits
        geneCopyNumbers_Accessible
        transcriptionUnitCopyNumbers_Accessible
        accessibleGenes	
        accessibleTranscriptionUnits
        ploidy
        
        linkingNumbers_minFreeEnergy %free energy mininum linking number
        supercoils                   %difference between linkingNumbers and free energy mininum linking number
        superhelicalDensity          %supercoils / free energy minimum linking number
        supercoiled                  %boolean (1 x 2) indicating whether or not each chromosome is properly supercoiled (eg. within some tolerance of the free energy minimum)
        
        damagedSites                 %integers (genome length x 4) indicating identities of bases which are damaged or which are adjacent to damaged bonds
        damagedSites_shifted_incm6AD %integers (genome length x 4) indicating identities of damaged bases/sites
        damagedSites_nonRedundant    %integers (genome length x 4) indicating identities of damaged bases/sites
        damagedSites_excm6AD         %integers (genome length x 4) indicating identities of damaged bases/sites
        gapSites3                    %boolean (genome length x 4) indicating positions of gap sites 3' to bases
        gapSites5                    %boolean (genome length x 4) indicating positions of gap sites 5' to bases
        abasicSites3                 %boolean (genome length x 4) indicating positions of abasic sites 3' to bases
        abasicSites5                 %boolean (genome length x 4) indicating positions of abasic sites 5' to bases
        damagedSugarPhosphates3      %integers (genome length x 4) indicating indices of damaged sugar-phosphate 3' to bases
        damagedSugarPhosphates5      %integers (genome length x 4) indicating indices of damaged sugar-phosphate 5' to bases
        damagedBases3                %integers (genome length x 4) indicating indices of damaged bases 3' to bases
        damagedBases5                %integers (genome length x 4) indicating indices of damaged bases 5' to bases
        strandBreaks3                %boolean (genome length x 4) indicating positions of strand breaks 3' to bases
        strandBreaks5                %boolean (genome length x 4) indicating positions of strand breaks 5' to bases
        intrastrandCrossLinks3       %boolean (genome length x 4) indicating positions of intrastrand cross links 3' to bases
        intrastrandCrossLinks5       %boolean (genome length x 4) indicating positions of intrastrand cross links 5' to bases
        hollidayJunctions3           %boolean (genome length x 4) indicating positions of holliday junctions 3' to bases
        hollidayJunctions5           %boolean (genome length x 4) indicating positions of holliday junctions 5' to bases
        
        singleStrandBreaks           %boolean (genome length x 4) indicating positions of single strand breaks (strand breaks excluding double strand breaks and strand breaks adjacent to gapSites)
        doubleStrandBreaks           %boolean (genome length x 4) indicating positions of double strand breaks        
        strandBreakClassification
        
        restrictableMunIRMSites      %boolean (genome length x 4) indicating positions of restrictable MunI R/M sites because they aren't methylated
        hemiunmethylatedMunIRMSites  %boolean (genome length x 4) indicating positions of hemi-unmethylated MunI R/M sites
        
        dryWeight                    %dry weight of this class' state properties
    end
    
    %references to objects
    properties
		genome		  %genomes
        compartment   %compartments
        gene          %genes
        
        metabolite    %metabolites
        transcript    %transcript
        rnaPolymerase %RNA polymerase
        
        dnaRepair     %DNA repair process
    end
    
    %constructor
    methods
        function this = Chromosome(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.CellState(wholeCellModelID, name);
        end
    end
    
    methods
        function storeObjectReferences(this, simulation)
			%this.genome = simulation.gene.genomes;%As M. genitalium only has one chromosome, while S. cerevisiae has 16 + 1 chromosomes
            this.compartment = simulation.compartment;
            this.gene = simulation.gene;
            this.metabolite = simulation.state('Metabolite');
            this.transcript = simulation.state('Transcript');
            this.rnaPolymerase = simulation.state('RNAPolymerase');
            this.dnaRepair = simulation.process('DNARepair');
        end
    end
    
    methods
        function initializeConstants(this, knowledgeBase, simulation)
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.CellState(knowledgeBase, simulation);
            
            %% import classes
            import edu.jiangnan.fmme.cell.sim.constant.ChromosomeSequence;
            
            %% sequence
			this.genome                 	= knowledgeBase.genomes;
			%this.sequence               	= {this.genome.sequence}';
			%this.sequenceLen            	= [this.genome.sequenceLength]';
			%this.sequenceGCContent      	= [this.genome.gcContent]';
			%this.ARSs                   	= {this.genome.genomeFeatures}';
			this.ARSs                   	= findobj(knowledgeBase.genomeFeatures, 'type', 'ARS');
			this.ARSStartpositions          = [this.ARSs.startCoordinate]';
			this.ARSTerpositions            = [this.ARSs.endCoordinate]';
			for i = 1:16
            this.sequence{i} = ChromosomeSequence(knowledgeBase.genomes(i).sequence);
            this.sequenceLen(i) = size(this.sequence{i}, 1);
			this.sequenceGCContent(i) = getGCContent(this.sequence{i});
			%Ensure that relaxedBasesPerTurn is defined such that the relaxed
            %linking number will be a whole number for this organism.
			this.relaxedBasesPerTurn(i)	 = length(this.sequence{i}) * (1 + this.equilibriumSuperhelicalDensity) / ...
				ceil(length(this.sequence{i}) / 10.5 * (1 + this.equilibriumSuperhelicalDensity));
			end
            %oriC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'oriC');
            %terC = findobj(knowledgeBase.genomeFeatures, 'wholeCellModelID', 'terC');
            %this.oriCPosition = oriC.startCoordinate;
            %this.terCPosition = terC.startCoordinate;
			
            %% transcription units
			for i = 1:16
			this.transcriptionUnitWholeCellModelIDs{i} = {knowledgeBase.genomes(i).transcriptionUnits.wholeCellModelID}';
			this.transcriptionUnitNames{i}  		   = {knowledgeBase.genomes(i).transcriptionUnits.name}';
			this.transcriptionUnitStartCoordinates{i}   = [knowledgeBase.genomes(i).transcriptionUnits.startCoordinate]';
			this.transcriptionUnitLengths{i}           = [knowledgeBase.genomes(i).transcriptionUnits.sequenceLength]';
			this.transcriptionUnitStrands{i}           = 2-[knowledgeBase.genomes(i).transcriptionUnits.direction]';
			end
            %% proteins 
			this.monomerIndexs_replisome     = zeros(0, 1);
			this.complexIndexs_replisome     = sort(simulation.state('ProteinComplex').getIndexs({
                'YKL045W_YIR008C_DIMER'
                'YOR217W_YJR068W_YNL290W_YOL094C_YBR087W_PENTAMER'
                'YKL017C_6MER'
                }));	
			
			%DNA ligase
			this.complexIndexs_ligase4	   = simulation.state('ProteinComplex').getIndexs({'DNA_LIGASE_IV'});%DNA ligase IV complex
			this.monomerIndexs_ligase1     = simulation.state('ProteinMonomer').getIndexs({'YDL164C_MONOMER'}); %DNA ligase 1
			%DNA Polymerase
			this.complexIndexs_dnaPolA	   = simulation.state('ProteinComplex').getIndexs({'DNA_POLYMERASE_A'});%alpha DNA polymerase:primase complex
			this.complexIndexs_dnaPolD	   = simulation.state('ProteinComplex').getIndexs({'DNA_POLYMERASE_D'});%delta DNA polymerase
			this.complexIndexs_dnaPolE	   = simulation.state('ProteinComplex').getIndexs({'DNA_POLYMERASE_E'});%Polymerase epsilon
			this.monomerIndexs_dnaPolG     = simulation.state('ProteinMonomer').getIndexs({'YOR330C_MONOMER'}); %Polymerase gamma
			this.complexIndexs_dnaPolZ	   = simulation.state('ProteinComplex').getIndexs({'DNA_POLYMERASE_Z'});%Polymerase zeta
			%G1/S DNA damage checkpoint 	
			this.monomerIndexs_BRE1        = simulation.state('ProteinMonomer').getIndexs({'YDL074C_MONOMER'});
			this.monomerIndexs_RAD9        = simulation.state('ProteinMonomer').getIndexs({'YDR217C_MONOMER'});
			this.monomerIndexs_RAD6        = simulation.state('ProteinMonomer').getIndexs({'YGL058W_MONOMER'});
			%G2/M DNA damage checkpoint
			this.monomerIndexs_CHK1        = simulation.state('ProteinMonomer').getIndexs({'YBR274W_MONOMER'});
            this.monomerIndexs_MUS81       = simulation.state('ProteinMonomer').getIndexs({'YDR386W_MONOMER'});
			this.monomerIndexs_DPB11       = simulation.state('ProteinMonomer').getIndexs({'YJL090C_MONOMER'});
			this.monomerIndexs_CHS5        = simulation.state('ProteinMonomer').getIndexs({'YLR330W_MONOMER'});
			%RNA Polymerase
			this.complexIndexs_rnaPol1	   = simulation.state('ProteinComplex').getIndexs({'RNA_POLYMERASE_I'});%DNA-directed RNA polymerase I complex
			this.complexIndexs_rnaPol2	   = simulation.state('ProteinComplex').getIndexs({'RNA_POLYMERASE_II'});%DNA-directed RNA polymerase II complex
			this.complexIndexs_rnaPol3	   = simulation.state('ProteinComplex').getIndexs({'RNA_POLYMERASE_III'});%DNA-directed RNA polymerase III complex
			
			
			
			this.monomerDNAFootprints = ceil([knowledgeBase.proteinMonomers.dnaFootprint]');
            this.complexDNAFootprints = ceil([knowledgeBase.proteinComplexs.dnaFootprint]');
            strandedValues = cell(3, 1);
            strandedValues{this.dnaStrandedness_ssDNA} = 'ssDNA';
            strandedValues{this.dnaStrandedness_dsDNA} = 'dsDNA';
            strandedValues{this.dnaStrandedness_xsDNA} = 'Either';
			this.complexDNAFootprints = ceil([knowledgeBase.proteinComplexs.dnaFootprint]');
			monomerStrandedness = {knowledgeBase.proteinMonomers.dnaFootprintBindingStrandedness}';
			complexStrandedness = {knowledgeBase.proteinComplexs.dnaFootprintBindingStrandedness}';
            monomerStrandedness(cellfun(@isempty, monomerStrandedness)) = {'dsDNA'};
            complexStrandedness(cellfun(@isempty, complexStrandedness)) = {'dsDNA'};
            [~, this.monomerDNAFootprintBindingStrandedness] = ismember(monomerStrandedness, strandedValues);
            [~, this.complexDNAFootprintBindingStrandedness] = ismember(complexStrandedness, strandedValues);
            monomerStrandedness = {knowledgeBase.proteinMonomers.dnaFootprintRegionStrandedness}';
            complexStrandedness = {knowledgeBase.proteinComplexs.dnaFootprintRegionStrandedness}';
            monomerStrandedness(cellfun(@isempty, monomerStrandedness)) = {'dsDNA'};
            complexStrandedness(cellfun(@isempty, complexStrandedness)) = {'dsDNA'};
            [~, this.monomerDNAFootprintRegionStrandedness] = ismember(monomerStrandedness, strandedValues);
            [~, this.complexDNAFootprintRegionStrandedness] = ismember(complexStrandedness, strandedValues);
            
            %% proteins which have ability to release other proteins bound to chromosomes
            state = findobj(knowledgeBase.states, 'wholeCellModelID', this.wholeCellModelID);
            reactions = state.reactions;
            reactionGlobalIndexs = [reactions.idx]';
            this.reactionWholeCellModelIDs = {reactions.wholeCellModelID}';
            this.reactionNames = {reactions.name}';
            
            this.reactionBoundMonomer = zeros(size(reactions));
            this.reactionBoundComplex = zeros(size(reactions));
            [i, j] = find(max(0, knowledgeBase.reactionProteinMonomerStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)));
            this.reactionBoundMonomer(j) = i;
            [i, j] = find(max(0, knowledgeBase.reactionProteinComplexStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)));
            this.reactionBoundComplex(j) = i;            
            this.reactionMonomerCatalysisMatrix = max(0, -knowledgeBase.reactionProteinMonomerStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)');
            this.reactionComplexCatalysisMatrix = max(0, -knowledgeBase.reactionProteinComplexStoichiometryMatrix(:, reactionGlobalIndexs, this.compartment.cytosolIndexs)');            
            this.reactionThresholds = ...
                sum(this.reactionMonomerCatalysisMatrix, 2) + ...
                sum(this.reactionComplexCatalysisMatrix, 2);
        end
    end
    
    methods
        %allocate memory
        function allocateMemory(this, numTimePoints)
            import edu.jiangnan.fmme.util.CircularSparseMat;
			for i = 1:16
			
            this.polymerizedRegions{i}            = CircularSparseMat([], [], [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.linkingNumbers{i}		          = CircularSparseMat([], [], [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);

            this.monomerBoundSites{i}             = CircularSparseMat([], [], [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.complexBoundSites{i}    	      = CircularSparseMat([], [], [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            
            this.gapSites{i}		              = CircularSparseMat([], false(0,1), [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.abasicSites{i}	                  = CircularSparseMat([], false(0,1), [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.damagedSugarPhosphates{i} 	      = CircularSparseMat([], [],         [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.damagedBases{i}                  = CircularSparseMat([], [],         [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.intrastrandCrossLinks{i}         = CircularSparseMat([], [],         [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.strandBreaks{i}                  = CircularSparseMat([], false(0,1), [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
            this.hollidayJunctions{i}             = CircularSparseMat([], false(0,1), [this.sequenceLen(i), this.nCompartments, numTimePoints], 1);
			end
			this.segregated	              	      = false(1, 1, numTimePoints);
            this.invalidate();
        end
    end
    
    methods
        function initialize(this)
            this.allocateMemory(1);
			for i = 1:16
            this.polymerizedRegions{i}(1, this.strandIndexs_ch1) = this.sequenceLen(i);
			%s	  											  = [this.sequenceLen] / [this.relaxedBasesPerTurn]' * (1 + this.equilibriumSuperhelicalDensity);
            %this.linkingNumbers(1, this.strandIndexs_ch1)     = {s(4,:)};
			this.linkingNumbers{i}(1, this.strandIndexs_ch1)	 = this.sequenceLen(i) / this.relaxedBasesPerTurn(i)' * (1 + this.equilibriumSuperhelicalDensity);
			end
            this.invalidate();
        end
    end

    %general methods which query, but do not modify state
    methods
        %- positionsStrands must be a 2 column vector
        %- lengths must have same number of rows as positionsStrands
        function [tfs, idxs, positionsStrands, lengths] = sampleAccessibleRegions(this, i, ...
                nSites, weights, positionsStrands, lengths, ...
                bindingMonomers, bindingComplexs, ...
                isPositionsStrandFootprintCentroid, ignoreDamageFilter, returnOverlappingRegions, returnExtentAccessible, ...
                checkRegionSupercoiled)
				
            if nargin < 13 %change from 12 to 13
                checkRegionSupercoiled = false;
            end
            
            if isempty(weights)
                weights = ones(size(positionsStrands, 1), 1);
            end
			
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness] = this.getDNAFootprint(i, bindingMonomers, bindingComplexs);
            
            idxs = [];
			
            while any(weights) && numel(idxs) < nSites
			
			
                %sample sites
                nMoreSites = min(max(2 * (nSites - numel(idxs)), 10), nnz(weights));
                selectedSites = this.randStream.randsample(numel(weights), nMoreSites, false, weights);
				
                weights(selectedSites) = 0;
                
                %determine if selected sites are accessible
                [tmpTfs, ~, ~, extents] = this.isRegionAccessible(i, positionsStrands(selectedSites, :), lengths(selectedSites), ...
                    bindingMonomers, bindingComplexs, isPositionsStrandFootprintCentroid, ignoreDamageFilter, ...
                    returnExtentAccessible, true, checkRegionSupercoiled);
				
				newIdxs = selectedSites(tmpTfs);
                lengths(selectedSites(~tmpTfs)) = 0;
                lengths(selectedSites(tmpTfs)) = extents;
                
                %if proteins binding, make sure proteins won't sterically overlap
                if ~returnOverlappingRegions
                    [idxs, newIdxs] = this.excludeOverlappingRegions(...
                        i, idxs, newIdxs, positionsStrands, lengths, ...
                        footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
                        footprintBindingStrandedness == this.dnaStrandedness_dsDNA);
                end
                
                if numel(newIdxs) > nSites - numel(idxs)
                    newIdxs = newIdxs(1:nSites - numel(idxs));
                end
                idxs = [idxs; newIdxs];
            end
            
            %format output
            idxs = sort(idxs);
            tfs = false(size(positionsStrands, 1), 1);
            tfs(idxs) = true;
            positionsStrands = positionsStrands(idxs, :);
            lengths = lengths(idxs, :);
        end
        %randomly select among accessible sites (which if sequence seq is
        %specified, have this sequence) with probability probOrNSites (or if
        %probOrNSites > 1, randomly select probOrNSites sites)
        function positionsStrands = sampleAccessibleSites(this, i, prob, nSites, seq)
		
		%for convenience
            dnaLength = this.sequenceLen(i);
            posStrnds = find(this.polymerizedRegions{i});
            nStrands = max(posStrnds(:, 2));
            seqLen = numel(seq);
            
            %estimate total number of accessible sites
            [~, boundMonomers] = find(this.monomerBoundSites{i});
            [~, boundComplexs] = find(this.complexBoundSites{i});
            
            nAccessibleSites  = ...
                collapse(this.polymerizedRegions{i} ) ...
                - sum(this.monomerDNAFootprints (boundMonomers , 1)) ...
                - sum(this.complexDNAFootprints (boundComplexs , 1)) ...
                - nnz(this.damagedSites{i});
            %calculate number of sites to select
            if ~isempty(prob)
                nGC = sum(seq == 'G' | seq == 'C');
				
                nSites = min(nSites, this.randStream.stochasticRound(...
                    nAccessibleSites * prob * ...
                    (this.sequenceGCContent(i)/2)^nGC * ((1-this.sequenceGCContent(i))/2)^(seqLen-nGC)));
            end
            if nSites == 0
                positionsStrands = zeros(0, 2);
                return;
            end
            
            %randomly select nSites undamaged sites with sequence equal to seq
            positionsStrands = zeros(0, 2);
            iter = 0;
            while iter < 15
                %iterations
                iter = iter + 1;
                
                %more sites
                nMoreSites = nSites - size(positionsStrands, 1);
                nMoreSites = max(2 * nMoreSites, nMoreSites + 10);
                
                %pick random positions and strand
                positions = ceil(dnaLength(i) * this.randStream.rand(nMoreSites, 1));
                strands   = ceil(nStrands * this.randStream.rand(nMoreSites, 1));
                
                %throw away positions that aren't inaccessible
                if isempty(seq)
                    [~, idxs] = this.isRegionAccessible(i, [positions strands], 1, [], [], true, [], false, false);
                else
                    dir = 2 * (mod(strands, 2) == 1) - 1;
                    pos = 0:seqLen(i) - 1;
                    
                    subsequences = this.sequence.subsequence(...
                        positions(:, ones(1, seqLen(i))) + ...
                        dir(:, ones(1, seqLen(i))) .* pos(ones(nMoreSites, 1), :), ...
                        strands);
                    
                    if size(seq, 1) == 1
                        if isscalar(seq)
                            idx = find(subsequences == seq);
                        elseif size(seq, 2) == 2
                            idx = find(subsequences(:, 1) == seq(1) & subsequences(:, 2) == seq(2));
                        else
                            idx = find(all(subsequences == seq(ones(size(subsequences, 1), 1), :), 2));
                        end
                    else
                        idx = find(ismember(subsequences, seq, 'rows'));
                    end
                    positions = positions(idx, :);
                    strands = strands(idx, :);
                    
                    [~, idxs] = this.isRegionAccessible(i, [positions strands], seqLen(i), [], [], true, [], false, false);
                end
                
                %append to list of valid positions and strands
                positionsStrands = edu.jiangnan.fmme.util.SparseMat.unique_subs(...
                    [positionsStrands; positions(idxs) strands(idxs)], [dnaLength(i) this.nCompartments]);
                
                %if more sites that requested, randomly select
                if size(positionsStrands, 1) >= nSites
                    positionsStrands = this.randStream.randomlySelectNRows(positionsStrands, nSites);
                    break;
                end
            end
            
            %sort
            if nSites > 1
                positionsStrands = edu.jiangnan.fmme.util.SparseMat.sort_subs(positionsStrands, [dnaLength(i) this.nCompartments]);
            end
        end
        
        %Filters a list of sites (positionsStrands -- strands and positions along
        %strands), returning only those sites which are accessible to the
        %queried protein monomers and complexes (that is site which have been
        %polymerized, which aren't damaged, and which are either not bound by
        %protein, or are bound by proteins which at least one of the query
        %protein can release from the chromosome). Query proteins are indicated
        %by their indices within simulation.matureMatureIndexs
        %(monomerIndexs) and simulation.matureIndexs (complexIndexs).
        %
        %If monomerIndexs and complexIndexs are null, or aren't defined, this
        %function only returns sites that have been polymerized, aren't damaged,
        %and aren't bound by any proteins.
        %
        %Also returns the indices (idxs) of the query sites which are
        %accessible, and a boolean vector (tfs) indicated whether or not each
        %query site is accessible
        %
        %Returns true/false if regions defined by positionsStrands and lengths
        %are accessible/inaccessible.
        %
        %true  ==> region is accessible
        %false ==> region is not accessible
        function [tfs, idxs, positionsStrands, lengths] = isRegionAccessible(this, i, ...
                positionsStrands, lengths, ...
                bindingMonomers, bindingComplexs, ...
                isPositionsStrandFootprintCentroid, ignoreDamageFilter, ...
                returnExtent, returnOverlappingRegions, checkRegionSupercoiled)

            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            %% default values of inputs
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput', 'Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if isempty(ignoreDamageFilter)
                ignoreDamageFilter = zeros(size(positionsStrands, 1), 0);
            elseif size(ignoreDamageFilter, 1) <= 1
                ignoreDamageFilter = ignoreDamageFilter(ones(size(positionsStrands, 1), 1), :);
            end
            if nargin < 10	%9
                returnOverlappingRegions = true;
            end
            if nargin < 11	%10
                checkRegionSupercoiled = false;
            end
            
            %% indices of complexs and monomers which at least one of the the
            %query proteins can release, and therefore shouldn't be cause for a
            %site to be filtered out if they are bound to a site
            [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(i, bindingMonomers, bindingComplexs);
            
			
            %% DNA footprint of binding proteins,
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness, footprintRegionStrandedness] = this.getDNAFootprint(i, bindingMonomers, bindingComplexs);
            
            %adding query positions if input positions are the centroid of
            %the binding proteins' DNA footpring
            queryPositionsStrands = positionsStrands;
            if isPositionsStrandFootprintCentroid
                strnd = queryPositionsStrands(:, 2);
                
                queryPositionsStrands(mod(strnd, 2) == 1 & lengths > 0, 1) = queryPositionsStrands(mod(strnd, 2) == 1 & lengths > 0, 1) - footprint5Prime;
                queryPositionsStrands(mod(strnd, 2) == 0 & lengths > 0, 1) = queryPositionsStrands(mod(strnd, 2) == 0 & lengths > 0, 1) - footprint3Prime;
                queryPositionsStrands(mod(strnd, 2) == 1 & lengths < 0, 1) = queryPositionsStrands(mod(strnd, 2) == 1 & lengths < 0, 1) + footprint3Prime;
                queryPositionsStrands(mod(strnd, 2) == 0 & lengths < 0, 1) = queryPositionsStrands(mod(strnd, 2) == 0 & lengths < 0, 1) + footprint5Prime;
                
                ignoreDamageFilter(mod(strnd, 2) == 1 & lengths > 0, :) = ignoreDamageFilter(mod(strnd, 2) == 1 & lengths > 0, :) + footprint5Prime;
                ignoreDamageFilter(mod(strnd, 2) == 0 & lengths > 0, :) = ignoreDamageFilter(mod(strnd, 2) == 0 & lengths > 0, :) + footprint3Prime;
                ignoreDamageFilter(mod(strnd, 2) == 1 & lengths < 0, :) = ignoreDamageFilter(mod(strnd, 2) == 1 & lengths < 0, :) - footprint3Prime;
                ignoreDamageFilter(mod(strnd, 2) == 0 & lengths < 0, :) = ignoreDamageFilter(mod(strnd, 2) == 0 & lengths < 0, :) - footprint5Prime;
            elseif lengths < 0
                queryPositionsStrands(:, 1) = queryPositionsStrands(:, 1) + footprint - 1;
            end
            queryLengths = lengths + sign(lengths) * (footprint - 1);
            
            %% filter query sites
            switch footprintRegionStrandedness
                case this.dnaStrandedness_dsDNA
                    [~, ~, ~, polymerized] = this.isRegionDoubleStranded(i, queryPositionsStrands, queryLengths, true, checkRegionSupercoiled);
                case this.dnaStrandedness_ssDNA
                    [~, ~, ~, polymerized] = this.isRegionSingleStranded(i, queryPositionsStrands, queryLengths, true);                
                case this.dnaStrandedness_xsDNA
                    %TODO: calculate extent
                    polymerized = queryLengths;
                otherwise
                    throw(MException('Chromosome:invalidInput','Invalid footprint strandedness'));
            end
			
            [~, ~, ~, proteinFree] = this.isRegionProteinFree(i, queryPositionsStrands, queryLengths, true, ...
                releasableMonomerIndexs, releasableComplexIndexs, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ...
                footprintRegionStrandedness  == this.dnaStrandedness_xsDNA);
	
            [~, ~, ~, undamaged] = this.isRegionUndamaged(queryPositionsStrands, queryLengths, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ignoreDamageFilter, true, ...
                footprintRegionStrandedness  == this.dnaStrandedness_xsDNA);
            
            extents = sign(lengths) .* max(0, (min(abs([proteinFree polymerized undamaged]), [], 2) - (footprint - 1)));
            
            %% Make sure regions don't sterically overlap
            if ~returnOverlappingRegions
                if returnExtent
                    tfs = find(extents ~= 0);
                else
                    tfs = find(lengths == extents);
                end
                [~, idxs] = this.excludeOverlappingRegions(i, [], tfs, positionsStrands, extents, ...
						footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
						footprintBindingStrandedness == this.dnaStrandedness_dsDNA);
					extents(setdiff(1:end, idxs)) = 0;
            end
             %- extract sites which pass filter
            %- reshape indices of sites which pass filter
            %- construct boolean indicating which sites pass filter
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = lengths == extents;
            end
            idxs = find(tfs);
            idxs = idxs(:);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
		end
    end
    
    %private methods for query state by region
    methods
        function [tfs, idxs, positionsStrands, lengths] = isRegionSingleStranded(this, i, positionsStrands, lengths, returnExtent)
			if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            [polPosStrnds, polLens] = find(this.polymerizedRegions{i});
            
            %only 1 chromosome
            if size(positionsStrands, 2) == 2 && all(polLens == this.sequenceLen(i)) && mod(numel(polLens), 2) == 0 && all(polPosStrnds(:, 2) == (1:numel(polLens))')
                if returnExtent
                    tfs = true(size(positionsStrands, 1), 1);
                    idxs = (1:numel(tfs))';
                    lengths(:) = 0;
                else
                    tfs = false(size(lengths));
                    idxs = zeros(0, 1);
                    positionsStrands = zeros(0, 2);
                    lengths = zeros(0, 1);
                end
                return;
            end
            %all other cases
            [tfs, idxs, positionsStrands, lengths] = this.isRegionPolymerized(i, positionsStrands, lengths, returnExtent, this.dnaStrandedness_ssDNA);
        end
        
        function [tfs, idxs, positionsStrands, lengths] = isRegionDoubleStranded(this, i, positionsStrands, lengths, returnExtent, checkRegionSupercoiled)
			if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if nargin < 6%5
                checkRegionSupercoiled = false;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            [polPosStrnds, polLens] = find(this.polymerizedRegions{i});
            if checkRegionSupercoiled
                supPosStrnds = find(this.supercoiled{i});
            end
            
            %only fully formed chromosomes
            if size(positionsStrands, 2) == 2 && all(polLens == this.sequenceLen(i)) && mod(numel(polLens), 2) == 0 && all(polPosStrnds(:, 2) == (1:numel(polLens))')
                if numel(polLens) == 2
                    if returnExtent
                        tfs = true(size(positionsStrands, 1), 1);
                        idxs = (1:numel(tfs))';
                        lengths(positionsStrands(:, 2) > 2) = 0;
                        if checkRegionSupercoiled
                            lengths(~ismembc(positionsStrands(:, 2), supPosStrnds(:, 2))) = 0;
                        end
                    else
                        tfs = positionsStrands(:, 2) <= 2;
                        if checkRegionSupercoiled
                            tfs = tfs & ismembc(positionsStrands(:, 2), supPosStrnds(:, 2));
                        end
                        idxs = find(tfs);
                        positionsStrands = positionsStrands(tfs, :);
                        lengths = lengths(tfs, 1);
                    end
                else
                    if returnExtent
                        tfs = true(size(positionsStrands, 1), 1);
                        idxs = (1:numel(tfs))';
                        if checkRegionSupercoiled
                            lengths(~ismembc(positionsStrands(:, 2), supPosStrnds(:, 2))) = 0;
                        end
                    else
                        tfs = true(size(positionsStrands, 1), 1);
                        if checkRegionSupercoiled
                            tfs = tfs & ismembc(positionsStrands(:, 2), supPosStrnds(:, 2));
                        end
                        idxs = find(tfs);
                        positionsStrands = positionsStrands(tfs, :);
                        lengths = lengths(tfs, 1);
                    end
                end
                return;
            end
            %all other cases
            [tfs, idxs, positionsStrands, lengths] = this.isRegionPolymerized(i, positionsStrands, lengths, returnExtent, this.dnaStrandedness_dsDNA, checkRegionSupercoiled);
		end
        
        function [tfs, idxs, positionsStrands, lengths] = isRegionPolymerized(this, i, ...
                positionsStrands, lengths, returnExtent, checkRegionStrandedness, checkRegionSupercoiled, polymerizedRegions)
			if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput','Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
				
            end
            if nargin < 8	%7
                polymerizedRegions = this.polymerizedRegions{i};
            end
            
            %if nargin >= 5 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_dsDNA
			if nargin >= 6 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_dsDNA
                [regions, regionLengths] = find(this.doubleStrandedRegions{i});
                %if nargin >= 6 && ~isempty(checkRegionSupercoiled) && checkRegionSupercoiled
				if nargin >= 7 && ~isempty(checkRegionSupercoiled) && checkRegionSupercoiled
                    tfs = this.supercoiled(regions);
                    regions = regions(tfs, :);
                    regionLengths = regionLengths(tfs, :);
                end
            %elseif nargin >= 5 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_ssDNA
             elseif nargin >= 6 && ~isempty(checkRegionStrandedness) && checkRegionStrandedness == this.dnaStrandedness_ssDNA   
				[regions, regionLengths] = find(this.singleStrandedRegions{i});
            else
                [regions, regionLengths] = find(polymerizedRegions);
            end
           
            %all strands either fully, or not at all synthesized
            if all(regionLengths == this.sequenceLen(i)) && size(regions, 2) == 2
                if size(regions, 1) == this.nCompartments
                    %all chromosomes fully synethesized
                    tfs = true(size(positionsStrands,1), 1);
                    idxs = (1:size(positionsStrands,1))';
                elseif isequal(regions(:,2), [1; 2])
                    %first chromosome fully synethesized
                    if returnExtent
                        tfs = true(size(positionsStrands,1), 1);
                        idxs = (1:size(positionsStrands,1))';
                        lengths(~ismembc(positionsStrands(:,2), regions(:, 2))) = 0;
                    else
                        tfs = ceil(positionsStrands(:,2)/2) == 1;
                        idxs = reshape(find(tfs), [], 1);
                        positionsStrands = positionsStrands(idxs, :);
                        lengths = lengths(idxs, :);
                    end
                else
                    %all strands either fully, or not at all synthesized
                    if returnExtent
                        tfs = true(size(positionsStrands,1), 1);
                        idxs = (1:size(positionsStrands,1))';
                        lengths(~ismembc(positionsStrands(:,2), regions(:, 2))) = 0;
                    else
                        tfs = ismembc(positionsStrands(:,2), regions(:, 2));
                        idxs = reshape(find(tfs), [], 1);
                        positionsStrands = positionsStrands(idxs, :);
                        lengths = lengths(idxs, :);
                    end
                end
                return;				
				%end
            end
            
            regionStartCoors = regions(:, 1);
            regionEndCoors = regionStartCoors  + regionLengths - 1;
            regionStrandTimes = [regions(:, 2:end) ones(size(regions, 1), size(positionsStrands,2) - size(regions, 2))];
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen(i)) +1;
            endCoors = startCoors + abs(lengths) - 1;
            strandTimes = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(regions, 2) - size(positionsStrands,2))];
            idxs = (1:size(positionsStrands, 1))';
            
            tmpIdxs = find(endCoors > this.sequenceLen(i));
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen(i))];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen(i)];
            endCoors(tmpIdxs) = this.sequenceLen(i);
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen(i)];
            endCoors = [endCoors; min(this.sequenceLen(i), endCoors(tmpIdxs) + this.sequenceLen(i))];
            startCoors(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands, 1), 1);
            if ndims(polymerizedRegions) == 2
                regionStrandTimeSiz = size(polymerizedRegions, 2);
            else
                regionStrandTimeSiz = [size(polymerizedRegions, 2) size(polymerizedRegions, 3)];
            end
            
            tfs = ...
                ((max(endCoors) >= regionStartCoors & max(endCoors) <= regionEndCoors) | ...
                (min(startCoors) >= regionStartCoors & min(startCoors) <= regionEndCoors) | ...
                (min(startCoors) < regionStartCoors & max(endCoors) > regionEndCoors)) & ...
                edu.jiangnan.fmme.util.SparseMat.ismember_subs(regionStrandTimes, strandTimes, regionStrandTimeSiz);
            regionStartCoors = regionStartCoors(tfs);
            regionEndCoors = regionEndCoors(tfs);
            regionStrandTimes = regionStrandTimes(tfs, :);
            
            if size(positionsStrands, 1) == numel(idxs)
                [~, idxs1, idxs2] = edu.jiangnan.fmme.util.SparseMat.unique_subs([startCoors endCoors], [this.sequenceLen(i); this.sequenceLen(i)]);
                for j = 1:numel(idxs1)
                    idxs3 = find(idxs2 == j);
                    startCoor = startCoors(idxs3(1));
                    endCoor = endCoors(idxs3(1));
                    len(i) = lengths(idxs(idxs3(1)));
                    if len(i) > 0
                        idxs3 = idxs3(idxs3 <= size(positionsStrands, 1) | extents(idxs(idxs3)) >= this.sequenceLen(i) - startCoor + 1);
                        tfs = startCoor >= regionStartCoors & startCoor <= regionEndCoors;
                        tmpregionEndCoors = regionEndCoors(tfs);
                        [tfs, idxs4] = edu.jiangnan.fmme.util.SparseMat.ismember_subs(strandTimes(idxs3, :), regionStrandTimes(tfs, :), regionStrandTimeSiz);
                        extents(idxs(idxs3(tfs))) = extents(idxs(idxs3(tfs))) + min(tmpregionEndCoors(idxs4(tfs)), endCoor) - startCoor + 1;
                    else
                        idxs3 = idxs3(idxs3 <= size(positionsStrands, 1) | -extents(idxs(idxs3)) >= endCoor);
                        tfs = endCoor >= regionStartCoors & endCoor <= regionEndCoors;
                        tmpregionStartCoors = regionStartCoors(tfs);
                        [tfs, idxs4] = edu.jiangnan.fmme.util.SparseMat.ismember_subs(strandTimes(idxs3, :), regionStrandTimes(tfs, :), regionStrandTimeSiz);
                        extents(idxs(idxs3(tfs))) = extents(idxs(idxs3(tfs))) - (endCoor - max(tmpregionStartCoors(idxs4(tfs)), startCoor) + 1);
                    end
                end
            else
                for j = 1:numel(startCoors)
                    if lengths(idxs(j)) > 0
                        if j > size(positionsStrands, 1) && extents(idxs(i)) < this.sequenceLen(i) - startCoors(idxs(j)) + 1
                            continue;
                        end
                        
                        tmpIdxs = find(...
                            startCoors(j) >= regionStartCoors & ...
                            startCoors(j) <= regionEndCoors & ...
                            ismember(regionStrandTimes, strandTimes(j, :), 'rows'), 1, 'first');
                        
                        if ~isempty(tmpIdxs)
                            extents(idxs(j)) = extents(idxs(j)) + min(regionEndCoors(tmpIdxs), endCoors(j)) - startCoors(j) + 1;
                        end
                    else
                        if j > size(positionsStrands, 1) && -extents(idxs(j)) < endCoors(idxs(j))
                            continue;
                        end
                        
                        tmpIdxs = find(...
                            endCoors(j) >= regionStartCoors & ...
                            endCoors(j) <= regionEndCoors & ...
                            ismember(regionStrandTimes, strandTimes(j, :), 'rows'), 1, 'first');
                        
                        if ~isempty(tmpIdxs)
                            extents(idxs(j)) = extents(idxs(j)) - (endCoors(j) - max(regionStartCoors(tmpIdxs), startCoors(j)) + 1);
                        end
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
		
        function [tfs, idxs, positionsStrands, lengths] = isRegionNotPolymerized(this, i, positionsStrands, lengths, returnExtent)
            if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            if ~all(lengths)
                throw(MException('Chromosome:invalidInput','Lengths must be non-zero integers'));
            elseif numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            
            unpolRgns = this.unpolymerizedRegions{i};
            if nnz(unpolRgns) == 0
                if returnExtent
                    tfs = true(size(positionsStrands,1), 1);
                    idxs = (1:size(positionsStrands,1))';
                    lengths =  zeros(size(positionsStrands,1), 1);
                else
                    tfs = false(size(positionsStrands,1), 1);
                    idxs = zeros(0, 1);
                    positionsStrands = zeros(0, ndims(positionsStrands));
                    lengths = zeros(0, 1);
                end
                return;
            end
            
            [regions, regionLengths] = find(unpolRgns);
            regionStartCoors = regions(:, 1);
            regionEndCoors = regionStartCoors  + regionLengths - 1;
            regionStrandTimes = [regions(: ,2:end) ones(size(regions, 1), size(positionsStrands,2) - size(regions, 2))];
            
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen(i)) +1;
            endCoors = startCoors + abs(lengths) - 1;
            strandTimes = [positionsStrands(: ,2:end) ones(size(positionsStrands, 1), size(regions, 2) - size(positionsStrands,2))];
            idxs = (1:size(positionsStrands, 1))';
            
            tmpIdxs = find(endCoors > this.sequenceLen(i));
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen(i))];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen(i)];
            endCoors(tmpIdxs) = this.sequenceLen(i);
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs, :)];
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen(i)];
            endCoors = [endCoors; min(this.sequenceLen(i), endCoors(tmpIdxs) + this.sequenceLen(i))];
            startCoors(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands, 1), 1);
            for j = 1:numel(startCoors)
                if lengths(idxs(j)) > 0
                    if j > size(positionsStrands, 1) && extents(idxs(j)) < this.sequenceLen(i) - startCoors(idxs(j)) + 1
                        continue;
                    end
                    
                    tmpIdxs = find(...
                        startCoors(j) >= regionStartCoors & ...
                        startCoors(j) <= regionEndCoors & ...
                        ismember(regionStrandTimes, strandTimes(j, :), 'rows'), 1, 'first');
                    
                    if ~isempty(tmpIdxs)
                        extents(idxs(j)) = extents(idxs(j)) + min(regionEndCoors(tmpIdxs), endCoors(j)) - startCoors(j) + 1;
                    end
                else
                    if j > size(positionsStrands, 1) && -extents(idxs(j)) < endCoors(idxs(j))
                        continue;
                    end
                    
                    tmpIdxs = find(...
                        endCoors(j) >= regionStartCoors & ...
                        endCoors(j) <= regionEndCoors & ...
                        ismember(regionStrandTimes, strandTimes(j, :), 'rows'), 1, 'first');
                    
                    if ~isempty(tmpIdxs)
                        extents(idxs(j)) = extents(idxs(j)) - (endCoors(j) - max(regionStartCoors(tmpIdxs), startCoors(j)) + 1);
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
		end

        %lengths must be non-negative integer
        %releasableMonomerIndexs and releasableComplexIndexs must be sorted
        function [tfs, idxs, positionsStrands, lengths, monomers, complexs] = isRegionProteinFree(this, i, ...
                positionsStrands, lengths, returnExtent, releasableMonomerIndexs, releasableComplexIndexs, bothStrands, allStrands)
			
			if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                monomers = zeros(0, 0);
                complexs = zeros(0, 0);
                return;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if nargin < 9	%8
                allStrands = false;
            end
			
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen(i)) + 1;
            endCoors = startCoors + abs(lengths) - 1;
            strnds = positionsStrands(:, 2);
            if bothStrands
                strnds = ceil(strnds/2);
            end
            if allStrands
                strnds(:) = 1;
            end
            idxs = (1:size(positionsStrands, 1))';
            tmpIdxs = find(endCoors > this.sequenceLen(i));
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen(i))];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen(i)];
			endCoors(tmpIdxs) = this.sequenceLen(i);
            strnds = [strnds; strnds(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen(i)];
            endCoors = [endCoors; min(this.sequenceLen(i), endCoors(tmpIdxs) + this.sequenceLen(i))];
            startCoors(tmpIdxs) = 1;
            strnds = [strnds; strnds(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs)];
            
            %initialize output
            extents  = abs(lengths);
            monomers = sparse(size(positionsStrands, 1), max(abs([lengths; 0])));
            complexs = sparse(size(positionsStrands, 1), max(abs([lengths; 0])));
            
            %bound monomers
            [subs, vals] = find(this.monomerBoundSites{i});
			
            if ~bothStrands
                dsTfs = this.monomerDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                subs = [
                    subs(~dsTfs, :)
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)-1
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)
                    ];
                vals = [
                    vals(~dsTfs)
                    vals(dsTfs)
                    vals(dsTfs)
                    ];
                [subs, order] = edu.jiangnan.fmme.util.SparseMat.sort_subs(subs, [this.sequenceLen(i) this.nCompartments]);
                vals = vals(order, :);
            end
            if allStrands
                tfs = ~ismembc(vals, this.monomerIndexs_replisome);
                subs = subs(tfs, :);
                vals = vals(tfs);
            end
            monomerStarts = subs(:, 1);
            monomerEnds   = monomerStarts + this.monomerDNAFootprints(vals, :) - 1;
            monomerStrnds = subs(:, 2);
            if bothStrands
                monomerStrnds = ceil(monomerStrnds/2);
            end
            if allStrands
                monomerStrnds(:) = 1;
            end
            
            tmpIdxs = find(monomerEnds > this.sequenceLen(i));
            monomerStarts = [monomerStarts; ones(numel(tmpIdxs), 1)];
            monomerEnds   = [monomerEnds;   monomerEnds(tmpIdxs, :) - this.sequenceLen(i)];
            monomerStrnds = [monomerStrnds; monomerStrnds(tmpIdxs, :)];
            vals = [vals; vals(tmpIdxs, :)];
            
            for j = 1:numel(startCoors)
                monMask = ...
                    (monomerStarts >= startCoors(j) & monomerStarts <= endCoors(j)) | ...
                    (monomerEnds >= startCoors(j) & monomerEnds <= endCoors(j)) | ...
                    (monomerStarts <= startCoors(j) & monomerEnds >= endCoors(j));
                tmpIdxs = find(monMask);
                if isempty(tmpIdxs); continue; end;
                
                tmpIdxs = tmpIdxs(monomerStrnds(monMask, 1) ==  strnds(j, 1));
                if isempty(tmpIdxs); continue; end;
                
                if lengths(idxs(j)) > 0
                    if j <= size(positionsStrands, 1)
                        monomerPos = monomerStarts(tmpIdxs) - startCoors(j) + 1;
                    else
                        monomerPos = monomerStarts(tmpIdxs) - startCoors(j) + 1 + (this.sequenceLen(i)-positionsStrands(idxs(j), 1));
                    end
                else
                    if j <= size(positionsStrands, 1)
                        monomerPos = endCoors(j) - monomerEnds(tmpIdxs) + 1;
                    else
                        monomerPos = endCoors(j) - monomerEnds(tmpIdxs) + 1 + positionsStrands(idxs(j), 1);
                    end
                end
                monomerPos = max(monomerPos, 1);
                
                extents(idxs(j)) = min([
                    extents(idxs(j));
                    monomerPos(~ismembc(vals(tmpIdxs), releasableMonomerIndexs))-1]);
                
                monomers(idxs(j), monomerPos) = vals(tmpIdxs); %#ok<SPRIX>
            end
            
            %bound complexes
            [subs, vals] = find(this.complexBoundSites{i});
            if ~bothStrands
                dsTfs = this.complexDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                subs = [
                    subs(~dsTfs, :)
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)-1
                    subs( dsTfs, 1) 2*ceil(subs(dsTfs, 2)/2)
                    ];
                vals = [
                    vals(~dsTfs)
                    vals( dsTfs)
                    vals( dsTfs)
                    ];
                [subs, order] = edu.jiangnan.fmme.util.SparseMat.sort_subs(subs, [this.sequenceLen(i) this.nCompartments]);
                vals = vals(order, :);
            end
            if allStrands
                tfs = ~ismembc(vals, this.complexIndexs_replisome);
                subs = subs(tfs, :);
                vals = vals(tfs);
            end
            complexStarts = subs(:, 1);
            complexEnds   = complexStarts + this.complexDNAFootprints(vals, :) - 1;
            complexStrnds = subs(:, 2);
            if bothStrands
                complexStrnds = ceil(complexStrnds/2);
            end
            if allStrands
                complexStrnds(:) = 1;
            end
            
            tmpIdxs = find(complexEnds > this.sequenceLen(i));
            complexStarts = [complexStarts; ones(numel(tmpIdxs), 1)];
            complexEnds   = [complexEnds;   complexEnds(tmpIdxs, :) - this.sequenceLen(i)];
            complexStrnds = [complexStrnds; complexStrnds(tmpIdxs, :)];
            vals = [vals; vals(tmpIdxs, :)];
            
            for j = 1:numel(startCoors)%change from i to j
                comMask = ...
                    (complexStarts >= startCoors(j) & complexStarts <= endCoors(j)) | ...
                    (complexEnds >= startCoors(j) & complexEnds <= endCoors(j)) | ...
                    (complexStarts <= startCoors(j) & complexEnds >= endCoors(j));
                tmpIdxs = find(comMask);
                if isempty(tmpIdxs); continue; end;
                
                tmpIdxs = tmpIdxs(complexStrnds(comMask, 1) ==  strnds(j, 1));
                if isempty(tmpIdxs); continue; end;
                
                if lengths(idxs(j)) > 0
                    if j <= size(positionsStrands, 1)
                        complexPos = complexStarts(tmpIdxs) - startCoors(j) + 1;
                    else
                        complexPos = complexStarts(tmpIdxs) - startCoors(j) + 1 + (this.sequenceLen(i)-positionsStrands(idxs(j), 1));
                    end
                else
                    if j <= size(positionsStrands, 1)
                        complexPos = endCoors(j) - complexEnds(tmpIdxs) + 1;
                    else
                        complexPos = endCoors(j) - complexEnds(tmpIdxs) + 1 + positionsStrands(idxs(j), 1);
                    end
                end
                complexPos = max(complexPos, 1);
                
                extents(idxs(j)) = min([
                    extents(idxs(j));
                    complexPos(~ismembc(vals(tmpIdxs), releasableComplexIndexs))-1]);
                
                complexs(idxs(j), complexPos) = vals(tmpIdxs); %#ok<SPRIX>
            end
            
            %format output
            extents = sign(lengths) .* extents;
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
        end
        %lengths must be non-negative integer
        function [tfs, idxs, positionsStrands, lengths] = isRegionUndamaged(this, ...
                positionsStrands, lengths, isEitherStrandDamaged, ignoreDamageFilter, returnExtent, isAnyStrandDamaged)
           %for i = 1:length(this.genome)
		   for i = 1:16
		   if isempty(positionsStrands)
                tfs = false(0, 1);
                idxs = zeros(0, 1);
                lengths = zeros(0, 1);
                return;
            end
            
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end
            if nargin < 7
                isAnyStrandDamaged = false;
            end
            
			[damagedPosStrnds, damages] = find(this.damagedSites{i});
            if isempty(damages)
                tfs = true(size(lengths));
                idxs = (1:numel(lengths))';
                return;
            end
            
            if isempty(ignoreDamageFilter)
                ignoreDamageFilter = zeros(size(positionsStrands, 1), 0);
            elseif size(ignoreDamageFilter, 1) <= 1
                ignoreDamageFilter = reshape(ignoreDamageFilter(ones(size(positionsStrands, 1), 1), :), size(positionsStrands, 1), []);
            end
                        
            damagedPositions   = damagedPosStrnds(:, 1);
            damagedStrandTimes = [damagedPosStrnds(:, 2:end) ones(size(damagedPosStrnds, 1), size(positionsStrands, 2) - size(damagedPosStrnds, 2))];
            strandTimes        = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(damagedPosStrnds, 2) - size(positionsStrands, 2))];
            
            if isAnyStrandDamaged
                strandTimes(:, 1) = 1;
                damagedStrandTimes(:, 1) = 1;
            elseif isEitherStrandDamaged
                strandTimes(:, 1) = ceil(strandTimes(:, 1) / 2);
                damagedStrandTimes(:, 1) = ceil(damagedStrandTimes(:, 1) / 2);
            end
                        
            startPositions = mod(positionsStrands(:, 1) + min(0, lengths+1) - 1, this.sequenceLen(i)) + 1;
            endPositions   = startPositions + abs(lengths) - 1;
            idxs = (1:size(positionsStrands,1))';
            
            tmpIdxs = find(endPositions > this.sequenceLen(i));
            strandTimes = [strandTimes; strandTimes(tmpIdxs)];
            startPositions = [startPositions; max(1, startPositions(tmpIdxs) - this.sequenceLen(i))];
            endPositions = [endPositions; endPositions(tmpIdxs) - this.sequenceLen(i)];
            endPositions(tmpIdxs) = this.sequenceLen(i);
            idxs = [idxs; idxs(tmpIdxs)];
            
            tmpIdxs = find(startPositions < 0);
            strandTimes = [strandTimes; strandTimes(tmpIdxs)];
            startPositions = [startPositions; startPositions(tmpIdxs) + this.sequenceLen(i)];
            endPositions = [endPositions; min(this.sequenceLen(i), endPositions(tmpIdxs) + this.sequenceLen(i))];
            startPositions(tmpIdxs) = 1;
            idxs = [idxs; idxs(tmpIdxs)];
            
            extents = zeros(size(positionsStrands,1), 1);
            for j = 1:numel(startPositions)%change from i to j
                tmpIdxs = ...
                    damagedPositions >= startPositions(j, 1) & ...
                    damagedPositions <= endPositions(j, 1) & ...
                    ismember(damagedStrandTimes, strandTimes(j, :), 'rows');
                
                if lengths(idxs(j)) > 0
                    if j > size(positionsStrands, 1) && extents(idxs(j)) < this.sequenceLen(i) - startPositions(idxs(j)) + 1
                        continue;
                    end
                    
                    if j <= size(positionsStrands, 1)
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(j), :);
                    else
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(j), :) - (this.sequenceLen(i) - startPositions(idxs(j)) + 1);
                    end
					
                    startPositions(j)
					damagedPositions
					tmpIgnoreDamageFilter
                    extent = min(setdiff(damagedPositions(tmpIdxs) - startPositions(j) + 1, tmpIgnoreDamageFilter)) - 1;
                    if isempty(extent)
                        extent = endPositions(j) - startPositions(j) + 1;
                    end
                    extents(idxs(j)) = extents(idxs(j)) + extent;
                else
                    if j <= size(positionsStrands, 1)
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(j), :);
                    else
                        tmpIgnoreDamageFilter = ignoreDamageFilter(idxs(j), :) + positionsStrands(idxs(j), 1);
                    end
                    
                    extent = min(setdiff(endPositions(j) - damagedPositions(tmpIdxs) + 1, tmpIgnoreDamageFilter - lengths(idxs(j)) - 1)) - 1;
                    if isempty(extent)
                        extents(idxs(j)) = extents(idxs(j)) -(endPositions(j) - startPositions(j) + 1);
                    else
                        extents(idxs(j)) = extents(idxs(j)) - extent;
                    end
                end
            end
            
            if returnExtent
                tfs = true(size(positionsStrands, 1), 1);
            else
                tfs = extents == lengths;
            end
            idxs = reshape(find(tfs), [], 1);
            positionsStrands = positionsStrands(idxs, :);
            lengths = extents(idxs, :);
			%end
        end
    end
    end
    %additional private methods used to query state
    methods
        %Converts from the base-pair centric view of the chromosomes that this
        %class uses to a strand-centric view where each column represents a
        %single strand.
        function varargout = getStrandView(this, i, outputs)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            if nargin < 3	%change from 2 to 3
                outputs = {
                    'polymerizedRegions';
                    'monomerBoundSites';
                    'complexBoundSites';
                    'gapSites';
                    'abasicSites';
                    'damagedSugarPhosphates';
                    'damagedBases';
                    'intrastrandCrossLinks';
                    'strandBreaks';
                    'hollidayJunctions';
                    };
            end
            if ~iscell(outputs)
                outputs = {outputs};
            end
            if numel(outputs) < nargout
                throw(MException('Chromosome:tooManyOutputs', 'Too many outputs requested'));
            end
            for k = 1:nargout
                tmp = this.(outputs{k}){i};
                switch outputs{k}
                    case 'polymerizedRegions'
						
                        [positionsStrands, lengths] = find(this.polymerizedRegions{i});
                        positionsStrands = [positionsStrands ones(size(positionsStrands,1), 3-size(positionsStrands,2))];
						
                        idx1 = find(positionsStrands(:,2) == this.strandIndexs_negative(1));
                        idx2 = find(positionsStrands(:,2) == this.strandIndexs_negative(2));
                        
                        starts1 = positionsStrands(idx1, 1);
                        starts2 = positionsStrands(idx2, 1);
						
                        ends1 = starts1 + lengths(idx1) - 1;
                        ends2 = starts2 + lengths(idx2) - 1;
                        times1 = positionsStrands(idx1, 3:end);
                        times2 = positionsStrands(idx2, 3:end);
                        
                        tmpStarts2 = zeros(0,1);
                        tmpEnds2 = zeros(0,1);
                        tmpTimes2 = zeros(0,size(times2,2));
						
                        for j = 1:numel(idx2)
                            idx = find(...
                                ((starts1 >= starts2(j) & starts1 <= ends2(j)) | ...
                                (ends1    >= starts2(j) & ends1   <= ends2(j)) | ...
                                (starts1  <= starts2(j) & ends1   >= ends2(j))) & ...
                                all(times1 == times2(j, :), 2));                                                                                   
							
                            tmpStarts2 = [tmpStarts2;
                                max(starts1(idx), starts2(j))];
                            tmpEnds2 = [tmpEnds2;
                                min(ends1(idx), ends2(j))];
                            tmpTimes2 = [tmpTimes2;
                                times1(idx, :)];
                            
                            starts1 = [starts1(setdiff(1:end, idx));
                                min([starts1(idx); starts2(j)])];
                            ends1 = [ends1(setdiff(1:end, idx));
                                max([ends1(idx); ends2(j)])];
                            times1 = [times1(setdiff(1:end, idx), :);
                                times2(j)];
                        end
                        tmpPosStrnds = [positionsStrands(setdiff(1:end, [idx1; idx2]), :);
                            starts1 this.strandIndexs_negative(ones(numel(starts1), 1)) times1;
                            tmpStarts2 this.strandIndexs_negative(2 * ones(numel(tmpStarts2), 1)) tmpTimes2];
                        tmpLengths = [lengths(setdiff(1:end, [idx1; idx2]), :);
                            ends1 - starts1 + 1;
                            tmpEnds2 - tmpStarts2 + 1;];            
						
                        tmp = CircularSparseMat(tmpPosStrnds, tmpLengths, size(tmp), 1);
                        tmp = this.mergeAdjacentRegions(tmp);
						
                    case 'linkingNumbers'
                        throw(MException('Chromosome:tooLazy', 'You''re going to have to write this yourself.'));
                    otherwise
                        [positionsStrands, lengths] = find(this.polymerizedRegions{i}(:, this.strandIndexs_negative(2), :));            
                        [tmpPosStrnds, tmpVals] = find(tmp);
						
                        positionsStrands = [positionsStrands ones(size(positionsStrands,1), 3-size(positionsStrands,2))];
                        tmpPosStrnds = [tmpPosStrnds ones(size(tmpPosStrnds,1), 3-size(tmpPosStrnds,2))];
					
                        for j = size(positionsStrands, 1)
                            idx1 = ...
                                tmpPosStrnds(:, 1) >= positionsStrands(j, 1) & ...
                                tmpPosStrnds(:, 1) <= positionsStrands(j, 1) + lengths(j)-1 & ...
                                tmpPosStrnds(:, 2) == this.strandIndexs_negative(1) & ...
                                all(tmpPosStrnds(:, 3:end) == positionsStrands(j, 3:end), 2);
                            idx2 = ...
                                tmpPosStrnds(:, 1) >= positionsStrands(j, 1) & ...
                                tmpPosStrnds(:, 1) <= positionsStrands(j, 1) + lengths(j)-1 & ...
                                tmpPosStrnds(:, 2) == this.strandIndexs_negative(2) & ...
                                all(tmpPosStrnds(:, 3:end) == positionsStrands{k}(j, 3:end), 2);
                            tmpPosStrnds(idx1, 2) = this.strandIndexs_negative(2);
                            tmpPosStrnds(idx2, 2) = this.strandIndexs_negative(1);
                        end
						
                        tmp = CircularSparseMat(tmpPosStrnds, tmpVals, size(tmp), 1);
					end
                varargout{k} = tmp;
            end            
        end        
        %Calculates the footprint of a group of proteins as the maximum of their
        %individual footprints. Also calculates the number of bases that the
        %footprint spans 3'- and 5'- to the centroid base of the footprint.
        function [footprint, footprint3Prime, footprint5Prime, bindingStrandedness, regionStrandedness] = getDNAFootprint(this, i, monomers, complexs)
			
			footprint = max([1;
                this.monomerDNAFootprints(monomers);
                this.complexDNAFootprints(complexs)]);
            
            [footprint3Prime, footprint5Prime] = this.calculateFootprintOverhangs(footprint);
            
            bindingStrandedness = [this.monomerDNAFootprintBindingStrandedness(monomers); this.complexDNAFootprintBindingStrandedness(complexs)];
            if isempty(bindingStrandedness)
                bindingStrandedness = this.dnaStrandedness_dsDNA;
            elseif ~isscalar(bindingStrandedness)
                if any(bindingStrandedness ~= bindingStrandedness(1))
                    throw(MException('Chromosome:invalidInput','Footprints must be of a single strandedness'));
                else
                    bindingStrandedness = bindingStrandedness(1);
                end
            end
            
            regionStrandedness = [this.monomerDNAFootprintRegionStrandedness(monomers); this.complexDNAFootprintRegionStrandedness(complexs)];
            if isempty(regionStrandedness)
                regionStrandedness = this.dnaStrandedness_dsDNA;
            elseif ~isscalar(regionStrandedness)
                if any(regionStrandedness ~= regionStrandedness(1))
                    throw(MException('Chromosome:invalidInput','Footprints must be of a single strandedness'));
                else
                    regionStrandedness = regionStrandedness(1);
                end
            end
        end
        %return sorted releasableMonomerIndexs, releasableComplexIndexs
        function [releasableMonomerIndexs, releasableComplexIndexs] = getReleasableProteins(this, i, bindingMonomers, bindingComplexs)
            releasableMonomerIndexs = this.reactionBoundMonomer .* (...
                sum(this.reactionMonomerCatalysisMatrix(:, bindingMonomers), 2) + ...
                sum(this.reactionComplexCatalysisMatrix(:, bindingComplexs), 2) >= ...
                this.reactionThresholds);
            
            if any(releasableMonomerIndexs)
                releasableMonomerIndexs = sort(releasableMonomerIndexs);
                if ~isempty(releasableMonomerIndexs)
                    releasableMonomerIndexs = releasableMonomerIndexs([diff(releasableMonomerIndexs) ~= 0; true], 1);
                end
            else
                releasableMonomerIndexs = 0;
            end
            
            releasableComplexIndexs = this.reactionBoundComplex .* (...
                sum(this.reactionMonomerCatalysisMatrix(:, bindingMonomers), 2) + ...
                sum(this.reactionComplexCatalysisMatrix(:, bindingComplexs), 2) >= ...
                this.reactionThresholds);
            if any(releasableComplexIndexs)
                releasableComplexIndexs = sort(releasableComplexIndexs);
                if ~isempty(releasableComplexIndexs)
                    releasableComplexIndexs = releasableComplexIndexs([diff(releasableComplexIndexs) ~= 0; true], 1);
                end
            else
                releasableComplexIndexs = 0;
            end
        end
        %Estimates regions that are accessible to a given protein species. That is
        %regions:
        %- that are polymerized, and have the strandedness that the protein requires
        %- that are not already bound by the specified protein species
        function [rgnPosStrnds, rgnLens] = getAccessibleRegions(this, i, monomerIdx, complexIdx, checkRegionSupercoiled)
			if nargin < 5 %change from 4 to 5
                checkRegionSupercoiled = false;
            end
            
            %find lengths and type of footprint
            [dnaFtpt, ~, ~, dnaFtptBindingStrandedness, dnaFtptRegionStrandedness] = this.getDNAFootprint(i, monomerIdx, complexIdx);
            
            %initialize excluded regions
            excPosStrnds = [];
            excLens = [];
            
            %exclude regions which are bound by proteins
            [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(i, monomerIdx, complexIdx);
            
            [monPosStrnds, mons] = find(this.monomerBoundSites{i});
            [cpxPosStrnds, cpxs] = find(this.complexBoundSites{i});
            
            iMonSS = find(~ismembc(mons, releasableMonomerIndexs) & this.monomerDNAFootprintBindingStrandedness(mons) == this.dnaStrandedness_ssDNA);
            iMonDS = find(~ismembc(mons, releasableMonomerIndexs) & this.monomerDNAFootprintBindingStrandedness(mons) == this.dnaStrandedness_dsDNA);
            iCpxSS = find(~ismembc(cpxs, releasableComplexIndexs) & this.complexDNAFootprintBindingStrandedness(cpxs) == this.dnaStrandedness_ssDNA);
            iCpxDS = find(~ismembc(cpxs, releasableComplexIndexs) & this.complexDNAFootprintBindingStrandedness(cpxs) == this.dnaStrandedness_dsDNA);
            
            excPosStrnds = [
                excPosStrnds
                monPosStrnds(iMonSS, :)
                monPosStrnds(iMonDS, 1) 2*ceil(monPosStrnds(iMonDS, 2)/2)-1
                monPosStrnds(iMonDS, 1) 2*ceil(monPosStrnds(iMonDS, 2)/2)
                cpxPosStrnds(iCpxSS, :)
                cpxPosStrnds(iCpxDS, 1) 2*ceil(cpxPosStrnds(iCpxDS, 2)/2)-1
                cpxPosStrnds(iCpxDS, 1) 2*ceil(cpxPosStrnds(iCpxDS, 2)/2)
                ];
            excLens = [
                excLens
                this.monomerDNAFootprints(mons(iMonSS))
                this.monomerDNAFootprints(mons(iMonDS))
                this.monomerDNAFootprints(mons(iMonDS))
                this.complexDNAFootprints(cpxs(iCpxSS))
                this.complexDNAFootprints(cpxs(iCpxDS))
                this.complexDNAFootprints(cpxs(iCpxDS))
                ];
            
            %exclude regions which are damaged
            dmgPosStrnds = find(this.damagedSites{i});
            
            excPosStrnds = [
                excPosStrnds;
                dmgPosStrnds];
            excLens = [
                excLens;
                ones(size(dmgPosStrnds, 1), 1)];
            
            %exclude other strand if proteins by both strands
            if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                excPosStrnds(:, 2) = ceil(excPosStrnds(:, 2) / 2);
            end
            
            %find polymerized regions
            switch dnaFtptRegionStrandedness
                case this.dnaStrandedness_ssDNA
                    if dnaFtptBindingStrandedness == this.dnaStrandedness_ssDNA
                        [polRgnPosStrnds, polRgnLens] = find(this.singleStrandedRegions{i});
                    else
                        throw(MException('Chromosome:error', 'unsupported strandedness: %s', dnaFtptRegionStrandedness));
                    end
                case this.dnaStrandedness_dsDNA
                    if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                        [polRgnPosStrnds, polRgnLens] = find(this.doubleStrandedRegions{i}(:, 1:2:end));
                        if checkRegionSupercoiled
                            tfs = this.supercoiled([polRgnPosStrnds(:, 1) 2 * polRgnPosStrnds(:, 2)]);
                            polRgnPosStrnds = polRgnPosStrnds(tfs, :);
                            polRgnLens = polRgnLens(tfs, :);
                        end
                    else
                        [polRgnPosStrnds, polRgnLens] = find(this.doubleStrandedRegions{i});
                        if checkRegionSupercoiled
                            tfs = this.supercoiled(polRgnPosStrnds);
                            polRgnPosStrnds = polRgnPosStrnds(tfs, :);
                            polRgnLens = polRgnLens(tfs, :);
                        end
                    end
                otherwise
                    throw(MException('Chromosome:error', 'unsupported strandedness: %s', dnaFtptRegionStrandedness));
            end
            
            %find difference of polymerized and excluded regions
            [rgnPosStrnds, rgnLens] = this.excludeRegions(i, polRgnPosStrnds, polRgnLens, excPosStrnds, excLens);
            
            %return strands
            if dnaFtptBindingStrandedness == this.dnaStrandedness_dsDNA
                rgnPosStrnds(:, 2) = 2 * rgnPosStrnds(:, 2) - 1;
            end
            
            %return only regions with length at least dnaFtpt
            idx = find(rgnLens >= dnaFtpt);
            rgnPosStrnds = rgnPosStrnds(idx, :);
            rgnLens = rgnLens(idx, :);
        end
        function value = getDamagedSites(this, i, includeBases, includeBonds, includeBase5Prime, includeBond5Prime, includeBase3Prime, includeBond3Prime, includeM6AD)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            if ~includeBases
                value = CircularSparseMat([], [], size(this.gapSites{i}), 1);
            else
                %damaged based
                if includeM6AD
                    value = this.damagedBases{i};
                else
                    [posStrnds, dmgs] = find(this.damagedBases{i});
                    posStrnds(dmgs == this.metabolite.m6ADIndexs, :) = [];
                    dmgs(dmgs == this.metabolite.m6ADIndexs, :) = [];
                    value = CircularSparseMat(posStrnds, dmgs, [this.sequenceLen(i) this.nCompartments], 1);
                end
                
                if nnz(this.gapSites{i})
                    value = value + this.gapSites{i};
                end
                if nnz(this.abasicSites{i})
                    value = value + this.abasicSites{i};
                end
                if nnz(this.damagedSugarPhosphates{i})
                    value = value + this.damagedSugarPhosphates{i};
                end
                if nnz(this.intrastrandCrossLinks{i})
                    value = value + this.intrastrandCrossLinks{i};
                    if includeBase5Prime
                        value = value + this.intrastrandCrossLinks5{i};
                    end
                    if includeBase3Prime
                        value = value + this.intrastrandCrossLinks3{i};
                    end
                end
            end
            
            if includeBonds
                if nnz(this.strandBreaks{i})
                    value = value + this.strandBreaks{i};
                    if includeBond5Prime
                        value = value + this.strandBreaks5{i};
                    end
                    if includeBond3Prime
                        value = value + this.strandBreaks3{i};
                    end
                end
                if nnz(this.hollidayJunctions{i})
                    value = value + this.hollidayJunctions{i};
                    if includeBond5Prime
                        value = value + this.hollidayJunctions5{i};
                    end
                    if includeBond3Prime
                        value = value + this.hollidayJunctions3{i};
                    end
                end
            end
        end
    end
    %additional methods which query, but do not modify state
    %unlike those above, these methods are very focused; they return information
    %about the Chromosome useful to only specific processes
    methods
        function [unmethylatedSites, hemimethylatedSites, methylatedSites, cleavedSites, inaccessibleRegions] = rmStatus(this, i, ...
                sites, methylatedPositions, restrictionPositions, bindingMonomers, bindingComplexs)
            import edu.jiangnan.fmme.util.SparseMat;
            
            warningState = warning('query', 'SparseMat:inefficient');
            warning('off', 'SparseMat:inefficient');
			for i = 1:16
            if nargin < 6	%change from 5 to 6
                bindingMonomers = [];
            end
            if nargin < 7	%change from 6 to 6
                bindingComplexs = [];
            end
            
            nonmethylatedPositions = [
                1:methylatedPositions(1)-1 methylatedPositions(1)+1:size(sites,2);
                1:methylatedPositions(2)-1 methylatedPositions(2)+1:size(sites,2)]';
            nonRestrictionPositions = [
                1:restrictionPositions(1)-1 restrictionPositions(1)+1:size(sites,2);
                1:restrictionPositions(2)-1 restrictionPositions(2)+1:size(sites,2)]';
            
            sitesChromosomesPositions = [
                size(sites, 1)
                size(sites, 2)
                size(this.damagedBases{i}, 2)/2
                size(this.damagedBases{i}, 3)
                ]';
            
            damagedRegions = ...
                (                this.damagedBases(                  sites(:,methylatedPositions(1)),                   1:2:end,:)~=this.metabolite.m6ADIndexs >0 & ...
                this.damagedBases(                  sites(:,methylatedPositions(1)),                   1:2:end,:)~=0                                       >0)   | ...
                (                this.damagedBases(                  sites(:,methylatedPositions(2)),                   2:2:end,:)~=this.metabolite.m6ADIndexs >0 & ...
                this.damagedBases(                  sites(:,methylatedPositions(2)),                   2:2:end,:)~=0                                       >0)   | ...
                xor(             this.strandBreaks(                  sites(:,restrictionPositions(1)),                  1:2:end,:), ...
                this.strandBreaks(                  sites(:,restrictionPositions(2)),                  2:2:end,:))                                               | ...
                collapse(reshape(this.damagedBases(          reshape(sites(:,nonmethylatedPositions(:,1)),        [],1), 1:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0   | ...
                collapse(reshape(this.damagedBases(          reshape(sites(:,nonmethylatedPositions(:,2)),        [],1), 2:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0   | ...
                this.damagedBases(                  sites(:,methylatedPositions(1)),                   2:2:end,:) ~=0                                            | ...
                this.damagedBases(                  sites(:,methylatedPositions(2)),                   1:2:end,:) ~=0                                            | ...
                collapse(reshape(this.strandBreaks(          reshape(sites(:,nonRestrictionPositions(1:end-1,1)),[],1), 1:2:end,:), sitesChromosomesPositions-[0 2 0 0]),2) >0    | ...
                collapse(reshape(this.strandBreaks(          reshape(sites(:,nonRestrictionPositions(2:end  ,2)),[],1), 2:2:end,:), sitesChromosomesPositions-[0 2 0 0]),2) >0    | ...
                this.strandBreaks(                  sites(:,restrictionPositions(1)),                  2:2:end,:) ~=0                                            | ...
                this.strandBreaks(                  sites(:,restrictionPositions(2)),                  1:2:end,:) ~=0                                            | ...
                collapse(reshape(this.intrastrandCrossLinks( reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.intrastrandCrossLinks( reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.damagedSugarPhosphates(reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.damagedSugarPhosphates(reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.abasicSites(           reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.abasicSites(           reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.gapSites(              reshape(sites,                                      [],1), 1:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.gapSites(              reshape(sites,                                      [],1), 2:2:end,:), sitesChromosomesPositions),2)           >0    | ...
                collapse(reshape(this.hollidayJunctions(     reshape(sites(:,1:end-1),                           [],1), 1:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0    | ...
                collapse(reshape(this.hollidayJunctions(     reshape(sites(:,2:end),                             [],1), 2:2:end,:), sitesChromosomesPositions-[0 1 0 0]),2) >0;
            
            methylatedSites = ...
                this.damagedBases{i}(sites(:, methylatedPositions(1)), 1:2:end,:) == this.metabolite.m6ADIndexs & ...
                this.damagedBases{i}(sites(:, methylatedPositions(2)), 2:2:end,:) == this.metabolite.m6ADIndexs;
            
            hemimethylatedSites = xor(...
                this.damagedBases{i}(sites(:, methylatedPositions(1)), 1:2:end,:) == this.metabolite.m6ADIndexs, ...
                this.damagedBases{i}(sites(:, methylatedPositions(2)), 2:2:end,:) == this.metabolite.m6ADIndexs);
            
            unmethylatedSites = ...
                this.damagedBases{i}(sites(:, methylatedPositions(1)), 1:2:end,:) ~= this.metabolite.m6ADIndexs & ...
                this.damagedBases{i}(sites(:, methylatedPositions(2)), 2:2:end,:) ~= this.metabolite.m6ADIndexs;
            
            cleavedSites = ...
                this.strandBreaks{i}(sites(:, restrictionPositions(1)), 1:2:end,:) & ...
                this.strandBreaks{i}(sites(:, restrictionPositions(2)), 2:2:end,:);
            
            siteLength = size(sites, 2);
            ignoreDamageFilter =  1:siteLength;
            sitesStrands = [
                sites(:, 1)   ones(size(sites, 1), 1)
                sites(:, 1) 3*ones(size(sites, 1), 1)];
            inaccessibleRegions = damagedRegions | ...
                (cleavedSites & ~unmethylatedSites) | ...
                ~reshape(this.isRegionAccessible(i, sitesStrands, siteLength, bindingMonomers, bindingComplexs, true, ignoreDamageFilter, false, true), [], this.nCompartments / 2);
                        
            methylatedSites     = ~inaccessibleRegions & methylatedSites;
            hemimethylatedSites = ~inaccessibleRegions & hemimethylatedSites;
            unmethylatedSites   = ~inaccessibleRegions & unmethylatedSites & ~cleavedSites;
            cleavedSites        = ~inaccessibleRegions & cleavedSites;            
            
            unmethylatedSites   = SparseMat(unmethylatedSites);
            hemimethylatedSites = SparseMat(hemimethylatedSites);
            methylatedSites     = SparseMat(methylatedSites);
            cleavedSites        = SparseMat(cleavedSites);
            inaccessibleRegions = SparseMat(inaccessibleRegions);
            
            if strcmp(warningState.state, 'on'); warning('on', 'SparseMat:inefficient'); end;
			end
		end
	end
    %methods which modify the state of this class, and possibly request
    %modifications to other parts of the simulation's state
    methods
        function sideEffects = setRegionUnwound(this, i, positions, lengths)
            sideEffects = edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect.empty(0, 1);
			
            if ~isequal(size(positions, 1), size(lengths, 1))
                throw(MException('Chromosome:invalidInput', 'positions and lengths must have same number of rows'));
            end
            if size(positions, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'positions must have 1 columns'));
            end
            if size(lengths, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'lengths must have 1 column'));
            end
            if ~any(lengths)
                return;
            end
            L(i) = this.sequenceLen(i);
            if any(positions > L(i))
                throw(MException('Chromosome:invalidInput', 'positions cannot wrap ORI'));
            end
            
            positions = positions(logical(lengths));
            lengths = lengths(logical(lengths));
            n = size(positions, 1);
            if ~all(this.isRegionDoubleStranded(i, [positions ones(n, 1)], lengths, false))
                throw(MException('Chromosome:invalidInput', 'regions must be double-stranded'));
            end
            if ~all(positions == 1 | positions == L(i) | ...
                    this.isRegionSingleStranded(i, [max(1, positions - 1) ones(n, 1)], 1, false) | ...
                    this.isRegionSingleStranded(i, [min(L(i), positions + 1) ones(n, 1)], 1, false))
                throw(MException('Chromosome:invalidInput',...
                    'unwinding must begin at either end of dsDNA or continue where it left off'));
            end
            if ~all(this.isRegionNotPolymerized(i, [
                    positions 3*ones(n, 1);
                    positions 4*ones(n, 1)],...
                    [lengths; lengths], false))
                throw(MException('Chromosome:invalidInput','chromosome 2 region cannot be polymerized'));
            end
            
            oldStrd = this.strandIndexs_ch1(2);
            newStrd = this.strandIndexs_ch2(2);
            
            directions = sign(lengths);
            positions = positions + min(0, lengths + 1);
            lengths = abs(lengths);
            
            for j = 1:n%change from i to j
                len = lengths(j);
                pos = positions(j,1);
                dir = directions(j);
                
                %if necessary, move region of initial negative strand of
                %chromosome 1 to chromosome 2
                this.monomerBoundSites{i}      = this.shiftStrandToNewChromosome(this.monomerBoundSites{i},      pos, len, oldStrd, newStrd);
                this.complexBoundSites{i}      = this.shiftStrandToNewChromosome(this.complexBoundSites{i},      pos, len, oldStrd, newStrd);
                this.gapSites{i}               = this.shiftStrandToNewChromosome(this.gapSites{i},               pos, len, oldStrd, newStrd);
                this.abasicSites{i}            = this.shiftStrandToNewChromosome(this.abasicSites{i},            pos, len, oldStrd, newStrd);
                this.damagedSugarPhosphates{i} = this.shiftStrandToNewChromosome(this.damagedSugarPhosphates{i}, pos, len, oldStrd, newStrd);
                this.damagedBases{i}           = this.shiftStrandToNewChromosome(this.damagedBases{i},           pos, len, oldStrd, newStrd);
                this.intrastrandCrossLinks{i}  = this.shiftStrandToNewChromosome(this.intrastrandCrossLinks{i},  pos, len, oldStrd, newStrd);
                this.strandBreaks{i}           = this.shiftStrandToNewChromosome(this.strandBreaks{i},           pos, len, oldStrd, newStrd);
                this.hollidayJunctions{i}      = this.shiftStrandToNewChromosome(this.hollidayJunctions{i},      pos, len, oldStrd, newStrd); 
                
                %update region of chromosomes that have been polymerized                
                [regionStartPositions{i}, regionLengths{i}] = find(this.polymerizedRegions{i}(:, oldStrd));
                regionStartPositions{i} = regionStartPositions{i}(:, 1);
                idx{i} = find(regionStartPositions{i} <= pos & regionStartPositions{i} + regionLengths{i} > pos);                
                this.polymerizedRegions{i}(regionStartPositions{i}(idx{i}), oldStrd) = pos - regionStartPositions{i}(idx{i});
                if dir == 1
                    if pos ~= regionStartPositions{i}(idx{i})
                        throw(MException('Chromosome:error', 'programmer error: unwinding bad region'));
                    end
                    this.linkingNumbers{i}(pos + len, 1:2) = this.linkingNumbers{i}(pos, oldStrd);
                    this.linkingNumbers{i}(pos, 1:2) = 0;
                end
                if (dir == 1 && len == regionLengths{i}(idx{i}) || ...
                    dir == -1 && pos == regionStartPositions{i}(idx{i})) && ...
                    abs(this.linkingNumbers{i}([pos oldStrd])) > 1e-6
                    throw(MException('Chromosome:invalidInput',...
                        'cannot completely unwind region with nonzero linking number'));
                end
                if regionLengths{i}(idx{i}) - (pos - regionStartPositions{i}(idx{i})) - len ~= 0
                    this.polymerizedRegions{i}(pos+len, oldStrd) = ...
                        regionLengths{i}(idx{i}) - (pos - regionStartPositions{i}(idx{i})) - len;
                end
                this.polymerizedRegions{i}(pos, newStrd) = len;
                
                this.mergeOwnAdjacentRegions(i);%need to be verified
            end
        end       
        function spmat = shiftStrandToNewChromosome(~, spmat, pos, len, oldStrd, newStrd)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            [subs, vals] = find(spmat);
            tfs = ...
                subs(:, 1) >= pos & ...
                subs(:, 1) <= pos + len - 1 & ...
                subs(:, 2) == oldStrd;
            subs(tfs, 2) = newStrd;
            spmat = CircularSparseMat(subs, vals, size(spmat), 1);
        end
        function sideEffects = setRegionPolymerized(this, i, positionsStrands, lengths)
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect;
            L(i) = this.sequenceLen(i);            
            if ~all(positionsStrands(:,2) == 1 | positionsStrands(:,2) == 2)
                throw(MException('Chromosome:invalidInput', 'positionsStrands must be a valid template strand (eg. 1 or 2)'));
            end
            if any(positionsStrands(:,1) + lengths -1 > L(i))
                throw(MException('Chromosome:invalidInput', 'positionsStrands cannot wrap ORI'));
            end
            if ~isequal(size(positionsStrands,1), size(lengths,1))
                throw(MException('Chromosome:invalidInput', 'positionsStrands and lengths must have same number of rows'));
            end
            if size(positionsStrands, 2) ~= 2
                throw(MException('Chromosome:invalidInput', 'positionsStrands must have 2 columns'));
            end
            if size(lengths, 2) ~= 1
                throw(MException('Chromosome:invalidInput', 'lengths must have 1 column'));
            end
            
            positionsStrands(:, 1) = positionsStrands(: ,1) + min(0, lengths + 1);
            lengths = abs(lengths);
            
            for j = 1:size(positionsStrands, 1)%change from i to j
                len(i) = lengths(j);
                pos = positionsStrands(j,1);
                tmpStrd = this.strandIndexs_template(positionsStrands(j,2), :);
                nonTmpStrd = this.strandIndexs_nonTemplate(positionsStrands(j,2), :);
                
                %if no polymerization requested, do nothing
                if len(i) == 0
                    continue;
                end

                %check that template exists
                if ~this.isRegionPolymerized(i, [pos tmpStrd], len(i), false)
                    throw(MException('Chromosome:error','cannot polymerize a region without a template'));
                end

                %check that strand hasn't been polymerized
                if ~this.isRegionNotPolymerized(i, [pos nonTmpStrd], len(i), false)
                    throw(MException('Chromosome:error','cannot polymerize a region that''s already been polymerized'));
                end

                %set region polymerized
                this.polymerizedRegions{i}([pos nonTmpStrd]) = len(i);
                this.linkingNumbers{i}([pos tmpStrd; pos nonTmpStrd]) = len(i) / this.relaxedBasesPerTurn(i);
                this.mergeOwnAdjacentRegions(i);
            end

            %side effects
            sideEffects = SimulationStateSideEffect.empty(0, 1);
        end
        %lengths must be non-negative integers
        function [releasedMonomers, releasedComplexs, sideEffects, tfs, idxs, positionsStrands, lengths] = ...
                setSiteProteinBound(this, i, positionsStrands, maxBindings, weights, binding_monomerIndexs, binding_complexIndexs, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs, ...
                isBindingStable, isPositionsStrandFootprintCentroid, lengths, isBindingProcessive, ignoreDamageFilter, checkRegionSupercoiled)
			
            if nargin < 15	%change from 14 to 15
                checkRegionSupercoiled = false;
            end
            
            if isBindingStable && numel(binding_monomerIndexs) + numel(binding_complexIndexs) > 1
                throw(MException('Chromosome:invalidInput', 'Can only bind 1 protein at a time'));
            end
            if isBindingStable && ~isBindingProcessive && ~all(lengths == 1)
                throw(MException('Chromosome:invalidInput','If any(lengths~=1) then if binding is stable it must be processive'));
            end
            if numel(lengths) == 1
                lengths = lengths(ones(size(positionsStrands, 1), 1), 1);
            end

            %randomly select among accessible sites
            [tfs, idxs, positionsStrands, lengths] = this.sampleAccessibleRegions(i, maxBindings, weights, positionsStrands, lengths, ...
                binding_monomerIndexs, binding_complexIndexs, isPositionsStrandFootprintCentroid, ignoreDamageFilter, ~isBindingStable, isBindingProcessive, checkRegionSupercoiled);
		
            %if positions are centroid, shift positions to start coordinate view
            [footprint, footprint3Prime, footprint5Prime, footprintBindingStrandedness, footprintRegionStrandedness] = this.getDNAFootprint(i, binding_monomerIndexs, binding_complexIndexs);
            
            releasePositionsStrands = positionsStrands;
            releasePositionsStrands(:, 1) = releasePositionsStrands(:, 1) + min(0, lengths + 1);
            if isPositionsStrandFootprintCentroid
                releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==1, 1) = releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==1, 1) - footprint5Prime;
                releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==0, 1) = releasePositionsStrands(mod(releasePositionsStrands(:,2),2)==0, 1) - footprint3Prime;
            end
            
            %remove proteins currently bound to selected accessible sites
            if footprintRegionStrandedness == this.dnaStrandedness_xsDNA
                [monomerPosStrands, monomers] = find(this.monomerBoundSites{i});
                [complexPosStrands, complexs] = find(this.complexBoundSites{i});
                [releasableMonomerIndexs, releasableComplexIndexs] = this.getReleasableProteins(i, binding_monomerIndexs, binding_complexIndexs);
                
                tmpIdxs = find(~ismembc(monomers, releasableMonomerIndexs));
                monomerPosStrands = monomerPosStrands(tmpIdxs, :);
                monomers = monomers(tmpIdxs, :);
                
                tmpIdxs = find(~ismembc(complexs, releasableComplexIndexs));
                complexPosStrands = complexPosStrands(tmpIdxs, :);
                complexs = complexs(tmpIdxs, :);
                
                [releasePositionsStrands, releaseLens] = this.excludeRegions(i, [
                    releasePositionsStrands(:, 1)     ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 2 * ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 3 * ones(size(releasePositionsStrands, 1), 1)
                    releasePositionsStrands(:, 1) 4 * ones(size(releasePositionsStrands, 1), 1)], ...
                    repmat(footprint + abs(lengths) - 1, 4, 1), ...
                    [monomerPosStrands; complexPosStrands], [this.monomerDNAFootprints(monomers); this.complexDNAFootprints(complexs)]);
            else
                releaseLens = footprint + abs(lengths) - 1;
            end
            [releasedMonomers, releasedComplexs, sideEffects] = ...
                this.setRegionProteinUnbound(i, releasePositionsStrands, releaseLens, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs, ...                
                footprintRegionStrandedness == this.dnaStrandedness_dsDNA, ...
                footprintBindingStrandedness == this.dnaStrandedness_dsDNA, ...
                false, false);
            
            %if binding is stable, bind protein to selected accessible sites
            if isBindingStable

                bindPositionsStrands = positionsStrands;
				
                if isBindingProcessive
                    bindPositionsStrands(:,1) = bindPositionsStrands(:,1) + lengths - sign(lengths);
                end
                if isPositionsStrandFootprintCentroid
                    bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==1, 1) = bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==1, 1) - footprint5Prime;
                    bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==0, 1) = bindPositionsStrands(mod(bindPositionsStrands(:,2),2)==0, 1) - footprint3Prime;
                end
                
                if ~isempty(binding_monomerIndexs)
                    this.monomerBoundSites{i}(bindPositionsStrands) = binding_monomerIndexs;
                    releasedMonomers(mainEffects_monomerIndexs == binding_monomerIndexs) = ...
                        releasedMonomers(mainEffects_monomerIndexs == binding_monomerIndexs) - ...
                        numel(idxs);
                else
                    this.complexBoundSites{i}(bindPositionsStrands) = binding_complexIndexs;
                    releasedComplexs(mainEffects_complexIndexs == binding_complexIndexs) = ...
                        releasedComplexs(mainEffects_complexIndexs == binding_complexIndexs) - ...
                        numel(idxs);
                end
            end
        end
        %change the identity of bound proteins
        function [releasedMonomers, releasedComplexs, sideEffects] = ...
                modifyBoundProtein(this, i, positionsStrands, newMonomers, newComplexs, ...
                mainEffects_monomerIndexs, mainEffects_complexIndexs)
             
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect;
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffectItem;
            import edu.jiangnan.fmme.util.countUnique;
			
            if any(all([newMonomers newComplexs], 2))
                throw(MException('Only 1 protein can be bound at each site'));
            end
            
            %get identities of currently bound proteins
            oldMonomers = this.monomerBoundSites{i}(positionsStrands);
            oldComplexs = this.complexBoundSites{i}(positionsStrands);
            
            %check that changing protein identify doesn't increase footprint
            %size (otherwise setSiteProteinBound should be used)
            oldFootprints = zeros(size(positionsStrands, 1));
            oldFootprints(oldMonomers ~= 0) = this.monomerDNAFootprints(oldMonomers(oldMonomers~=0));
            oldFootprints(oldComplexs ~= 0) = this.complexDNAFootprints(oldComplexs(oldComplexs~=0));
            
            newFootprints = zeros(size(positionsStrands, 1));
            newFootprints(newMonomers ~= 0) = this.monomerDNAFootprints(newMonomers(newMonomers~=0));
            newFootprints(newComplexs ~= 0) = this.complexDNAFootprints(newComplexs(newComplexs~=0));
            
            if any(newFootprints > oldFootprints)
                throw(MException('New protein footprints cannot be larger that old ones'));
            end
            
            %modify identities of bound proteins
            this.monomerBoundSites{i}(positionsStrands(newMonomers~=0, :)) = newMonomers(newMonomers~=0, :);
            this.complexBoundSites{i}(positionsStrands(newComplexs~=0, :)) = newComplexs(newComplexs~=0, :);
			
            %summarize bound/released monomers/complexs
            releasedMonomers = zeros(size(mainEffects_monomerIndexs));
            releasedComplexs = zeros(size(mainEffects_complexIndexs));
            
            oldMonomers = oldMonomers(oldMonomers~=0);
            [tfs, idxs] = ismember(oldMonomers, mainEffects_monomerIndexs);
			
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
			
			[idxs, counts] = countUnique(idxs);
            releasedMonomers(idxs) = releasedMonomers(idxs) + counts;
            
            newMonomers = newMonomers(newMonomers~=0);
            [tfs, idxs] = ismember(newMonomers, mainEffects_monomerIndexs);
			
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedMonomers(idxs) = releasedMonomers(idxs) - counts;
            
            oldComplexs = oldComplexs(oldComplexs~=0);
            [tfs, idxs] = ismember(oldComplexs, mainEffects_complexIndexs);
			
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedComplexs(idxs) = releasedComplexs(idxs) + counts;
            
            newComplexs = newComplexs(newComplexs~=0);
            [tfs, idxs] = ismember(newComplexs, mainEffects_complexIndexs);
			
            if ~all(tfs)
                throw(MException('Chromosome:invalidInput','All binding/releasing proteins must be mainEffects'));
            end
            [idxs, counts] = countUnique(idxs);
            releasedComplexs(idxs) = releasedComplexs(idxs) - counts;
            
            sideEffects = SimulationStateSideEffect.empty(0, 1);
        end
        function [releasedMonomers, releasedComplexs, sideEffects] = setRegionProteinUnbound(this, i, ...
                positionsStrands, lengths, mainEffects_monomerIndexs, mainEffects_complexIndexs, ...
                regionBothStrands, bindingBothStrands, suspendExternalStateUpdating, proteinIsDegraded)
				
			import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect;
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffectItem;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            import edu.jiangnan.fmme.util.countUnique;
            startCoors = mod(positionsStrands(:, 1) + min(0, lengths +1) - 1, this.sequenceLen(i)) + 1;
            endCoors = startCoors + abs(lengths) - 1;
			
            strndTimes = [positionsStrands(:, 2:end) ones(size(positionsStrands, 1), size(this.monomerBoundSites{i}, 3) - size(positionsStrands,2))];
            if bindingBothStrands
                strndTimes(:, 1) = ceil(strndTimes(:, 1) / 2); 
            end
            
            tmpIdxs = find(endCoors > this.sequenceLen(i));
            startCoors = [startCoors; max(1, startCoors(tmpIdxs) - this.sequenceLen(i))];
            endCoors = [endCoors; endCoors(tmpIdxs) - this.sequenceLen(i)];
            endCoors(tmpIdxs) = this.sequenceLen(i);
            strndTimes = [strndTimes; strndTimes(tmpIdxs, :)];
            
            tmpIdxs = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmpIdxs) + this.sequenceLen(i)];
            endCoors = [endCoors; min(this.sequenceLen(i), endCoors(tmpIdxs) + this.sequenceLen(i))];
            startCoors(tmpIdxs) = 1;
            strndTimes = [strndTimes; strndTimes(tmpIdxs, :)];
            
            %bound monomers
            [subs, vals] = find(this.monomerBoundSites{i});
            unbindingMonomers = false(size(vals));
            monomerStarts = subs(:, 1);
            monomerEnds   = monomerStarts + this.monomerDNAFootprints(vals, :) - 1;
            monomerStrndTimes = [subs(:, 2:end) ones(size(subs,1), size(strndTimes,2)-size(subs,2)+1)];
            idxs = (1:numel(vals))';
            if bindingBothStrands
                monomerStrndTimes(:, 1) = ceil(monomerStrndTimes(:, 1) / 2);
            elseif regionBothStrands
                tfs = this.monomerDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                idxs = [idxs(~tfs); idxs(tfs); idxs(tfs)];
                monomerStarts = [monomerStarts(~tfs); monomerStarts(tfs); monomerStarts(tfs)];
                monomerEnds = [monomerEnds(~tfs); monomerEnds(tfs); monomerEnds(tfs)];
                monomerStrndTimes = [
                    monomerStrndTimes(~tfs, :)
                    2*ceil(monomerStrndTimes(tfs, 1)/2)-1 monomerStrndTimes(tfs, 2:end)
                    2*ceil(monomerStrndTimes(tfs, 1)/2) monomerStrndTimes(tfs, 2:end)
                    ];
            end
            
            tmpIdxs = find(monomerEnds > this.sequenceLen(i));
            monomerStarts = [monomerStarts; ones(numel(tmpIdxs), 1)];
            monomerEnds   = [monomerEnds;   monomerEnds(tmpIdxs, :) - this.sequenceLen(i)];
            monomerStrndTimes = [monomerStrndTimes; monomerStrndTimes(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs, :)];
            
            for j = 1:numel(startCoors) %change from i to j 
                monMask = ...
                    (monomerStarts >= startCoors(j) & monomerStarts <= endCoors(j)) | ...
                    (monomerEnds >= startCoors(j) & monomerEnds <= endCoors(j)) | ...
                    (monomerStarts <= startCoors(j) & monomerEnds >= endCoors(j));
                tmpIdxs = find(monMask);
                if isempty(tmpIdxs); continue; end;
                
                if isvector(strndTimes)
                    tmpIdxs = tmpIdxs(monomerStrndTimes(monMask, 1) ==  strndTimes(j, 1));
                elseif isscalar(tmpIdxs)
                    tmpIdxs = tmpIdxs(all(monomerStrndTimes(monMask, :) == strndTimes(j, :)));
                else
                    tmpIdxs = tmpIdxs(ismember(monomerStrndTimes(monMask, :), strndTimes(j, :), 'rows'));
                end
                if isempty(tmpIdxs); continue; end;
                
                unbindingMonomers(idxs(tmpIdxs)) = true;
            end
            
            releasedMonomers = zeros(numel(mainEffects_monomerIndexs), 1);
            if isempty(unbindingMonomers) || ~any(unbindingMonomers)
                gblIdxs = [];
                counts = [];
            elseif isscalar(unbindingMonomers)
                gblIdxs = vals(unbindingMonomers);
                lclIdxs = find(gblIdxs == mainEffects_monomerIndexs);
                if isempty(lclIdxs)
                    counts = 1;
                else
                    releasedMonomers(lclIdxs, 1) = 1;
                    counts = zeros(0, 1);
                    gblIdxs = zeros(0, 1);
                end
            else
                [gblIdxs, counts] = countUnique(vals(unbindingMonomers, 1));
                [tfs, lclIdxs] = ismember(gblIdxs, mainEffects_monomerIndexs);
                
                releasedMonomers(lclIdxs(tfs)) = counts(tfs);
                
                tfs = tfs | gblIdxs == 0;
                gblIdxs = gblIdxs(~tfs);
                counts = counts(~tfs);
            end

            if isempty(gblIdxs)
                sideEffects_releasedMonomers = SimulationStateSideEffect.empty(0, 1);
            else
                sideEffects_releasedMonomers = SimulationStateSideEffect.empty(numel(gblIdxs), 0);
            end
            for j = 1:numel(gblIdxs) %change from i to j
                sideEffects_releasedMonomers(j, 1) = SimulationStateSideEffect([...
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', 'matureIndexs', gblIdxs(j), this.compartment.cytosolIndexs,  counts(j)); ...
                    SimulationStateSideEffectItem('ProteinMonomer', 'counts', 'boundIndexs',  gblIdxs(j), this.compartment.cytosolIndexs, -counts(j))]);
            end
            unbindingMonomerSubs = subs(unbindingMonomers, :);
            unbindingMonomerVals = vals(unbindingMonomers, :);
            this.monomerBoundSites{i} = CircularSparseMat(subs(~unbindingMonomers, :), vals(~unbindingMonomers, 1), [this.sequenceLen(i) this.nCompartments], 1);
            
            %bound complexes
            [subs, vals] = find(this.complexBoundSites{i});
            unbindingComplexs = false(size(vals));
            complexStarts = subs(:, 1);
            complexEnds   = complexStarts + this.complexDNAFootprints(vals, :) - 1;
            complexStrndTimes = [subs(:, 2:end) ones(size(subs,1), size(strndTimes,2)-size(subs,2)+1)];
            idxs = (1:numel(vals))';
            if bindingBothStrands
                complexStrndTimes(:, 1) = ceil(complexStrndTimes(:, 1) / 2);
            elseif regionBothStrands
                tfs = this.complexDNAFootprintBindingStrandedness(vals) == this.dnaStrandedness_dsDNA;
                idxs = [idxs(~tfs); idxs(tfs); idxs(tfs)];
                complexStarts = [complexStarts(~tfs); complexStarts(tfs); complexStarts(tfs)];
                complexEnds = [complexEnds(~tfs); complexEnds(tfs); complexEnds(tfs)];
                complexStrndTimes = [
                    complexStrndTimes(~tfs, :)
                    2*ceil(complexStrndTimes(tfs, 1)/2)-1 complexStrndTimes(tfs, 2:end)
                    2*ceil(complexStrndTimes(tfs, 1)/2) complexStrndTimes(tfs, 2:end)
                    ];
            end
            
            tmpIdxs = find(complexEnds > this.sequenceLen(i));
            complexStarts = [complexStarts; ones(numel(tmpIdxs), 1)];
            complexEnds   = [complexEnds;   complexEnds(tmpIdxs, :) - this.sequenceLen(i)];
            complexStrndTimes    = [complexStrndTimes; complexStrndTimes(tmpIdxs, :)];
            idxs = [idxs; idxs(tmpIdxs, :)];
            
            for j = 1:numel(startCoors)%change from i to j
                comMask = ...
                    (complexStarts >= startCoors(j) & complexStarts <= endCoors(j)) | ...
                    (complexEnds >= startCoors(j) & complexEnds <= endCoors(j)) | ...
                    (complexStarts <= startCoors(j) & complexEnds >= endCoors(j));
                tmpIdxs = find(comMask);
                if isempty(tmpIdxs); continue; end;
                
                if isvector(strndTimes)
                    tmpIdxs = tmpIdxs(complexStrndTimes(comMask, 1) ==  strndTimes(j, 1));
                elseif isscalar(tmpIdxs)
                    tmpIdxs = tmpIdxs(all(complexStrndTimes(comMask, :) == strndTimes(j, :)));
                else
                    tmpIdxs = tmpIdxs(ismember(complexStrndTimes(comMask, :), strndTimes(j, :), 'rows'));
                end
                if isempty(tmpIdxs); continue; end;
                
                unbindingComplexs(idxs(tmpIdxs)) = true;
            end
            
            releasedComplexs = zeros(numel(mainEffects_complexIndexs), 1);
            if isempty(unbindingComplexs) || ~any(unbindingComplexs)
                gblIdxs = [];
                counts = [];
            elseif isscalar(unbindingComplexs)
                gblIdxs = vals(unbindingComplexs);
                lclIdxs = find(gblIdxs == mainEffects_complexIndexs);
                if isempty(lclIdxs)
                    counts = 1;
                else
                    releasedComplexs(lclIdxs, 1) = 1;
                    counts = zeros(0, 1);
                    gblIdxs = zeros(0, 1);
                end
            else
                [gblIdxs, counts] = countUnique(vals(unbindingComplexs, 1));
                [tfs, lclIdxs] = ismember(gblIdxs, mainEffects_complexIndexs);
                
                releasedComplexs(lclIdxs(tfs)) = counts(tfs);
                tfs = tfs | gblIdxs == 0;
                gblIdxs = gblIdxs(~tfs);
                counts = counts(~tfs);
            end
            
            if isempty(gblIdxs)
                sideEffects_releasedComplexs = SimulationStateSideEffect.empty(0, 1);
            else
                sideEffects_releasedComplexs = SimulationStateSideEffect.empty(numel(gblIdxs), 0);
            end
            for j = 1:numel(gblIdxs)
                sideEffects_releasedComplexs(j, 1) = SimulationStateSideEffect([...
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', 'matureIndexs', gblIdxs(j), this.compartment.cytosolIndexs,  counts(j)); ...
                    SimulationStateSideEffectItem('ProteinComplex', 'counts', 'boundIndexs',  gblIdxs(j), this.compartment.cytosolIndexs, -counts(j))]);
            end
			
            unbindingComplexSubs = subs(unbindingComplexs, :);
            unbindingComplexVals = vals(unbindingComplexs, :);
            this.complexBoundSites{i} = CircularSparseMat(subs(~unbindingComplexs, :), vals(~unbindingComplexs, 1), [this.sequenceLen(i) this.nCompartments], 1);
            
            %effects on other states
            if ~suspendExternalStateUpdating
                this.updateExternalState(i, unbindingMonomerSubs, unbindingMonomerVals, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded);
            end
            
            %side effects
            sideEffects = [
                sideEffects_releasedMonomers;
                sideEffects_releasedComplexs];
        end
        function [positionsStrands, sideEffects] = stochasticallySetProteinUnbound(this, i, monomerIndex, complexIndex, ...
                rate, protectedPositionsStrands, protectedLengths, suspendExternalStateUpdating, proteinIsDegraded)
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            sideEffects = SimulationStateSideEffect.empty(0, 1);
            if ~isempty(monomerIndex)
                if ~isempty(complexIndex)
                    throw(MException('Chromosome:invalidInput', 'Can only unbind one protein at a time'));
                end
                [positionsStrands, proteins] = find(this.monomerBoundSites{i});
                unbinding = proteins == monomerIndex;
            elseif ~isempty(complexIndex)
                [positionsStrands, proteins] = find(this.complexBoundSites{i});
                unbinding = proteins == complexIndex;
            else
                positionsStrands = zeros(0, 2);
                return;
            end
            if ~any(unbinding)
                positionsStrands = zeros(0, 2);
                return;
            end
            
            %randomly select bound protein to release at specified rate
            if isfinite(rate)
                unbinding(unbinding) = this.randStream.rand(sum(unbinding), 1) < rate;
                if ~any(unbinding)
                    positionsStrands = zeros(0, 2);
                    return;
                end
            end
            
            %exclude protected sites
            if ~isempty(protectedPositionsStrands)
                footprint = this.getDNAFootprint(i, monomerIndex, complexIndex);
                for j = 1:size(protectedPositionsStrands, 1) %change fro i to j
                    unbndPosStrnds = positionsStrands(unbinding, :);
                    
                    tfs = ...
                        (unbndPosStrnds(:, 1)                 >= protectedPositionsStrands(j, 1) & ...
                         unbndPosStrnds(:, 1)                 <= protectedPositionsStrands(j, 1) + protectedLengths(j) - 1) | ...
                        (unbndPosStrnds(: ,1) + footprint - 1 >= protectedPositionsStrands(j, 1) & ...
                         unbndPosStrnds(:, 1) + footprint - 1 <= protectedPositionsStrands(j, 1) + protectedLengths(j) - 1) | ...
                        (unbndPosStrnds(:, 1)                 <= protectedPositionsStrands(j, 1) & ...
                         unbndPosStrnds(:, 1) + footprint - 1 >= protectedPositionsStrands(j, 1) + protectedLengths(j) - 1);
                    
                    tfs(tfs) = unbndPosStrnds(tfs, 2) == protectedPositionsStrands(j, 2);
                    
                    unbinding(unbinding) = ~tfs;
                    
                    if all(tfs)
                        positionsStrands = zeros(0, 2);
                        return;
                    end
                end
            end
            
            %update bound proteins
            if ~isempty(monomerIndex)
                unbindingMonomerSubs = positionsStrands(unbinding, :);
                unbindingMonomerVals = proteins(unbinding, 1);
                unbindingComplexSubs = zeros(0, 2);
                unbindingComplexVals = zeros(0, 1);
                this.monomerBoundSites{i} = CircularSparseMat(positionsStrands(~unbinding, :), proteins(~unbinding, 1), [this.sequenceLen(i) this.nCompartments], 1);
            else
                unbindingMonomerSubs = zeros(0, 2);
                unbindingMonomerVals = zeros(0, 1);
                unbindingComplexSubs = positionsStrands(unbinding, :);
                unbindingComplexVals = proteins(unbinding, 1);
                this.complexBoundSites{i} = CircularSparseMat(positionsStrands(~unbinding, :), proteins(~unbinding, 1), [this.sequenceLen(i) this.nCompartments], 1);
            end
            positionsStrands = positionsStrands(unbinding, :);
            
            %effects on other states
            if ~suspendExternalStateUpdating
                this.updateExternalState(i, unbindingMonomerSubs, unbindingMonomerVals, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded);
            end
        end
        function updateExternalState(this, i, ~, ~, unbindingComplexSubs, unbindingComplexVals, proteinIsDegraded)
            %posStrnds = unbindingComplexSubs(ismembc(unbindingComplexVals, this.complexIndexs_rnaPolymerase), :);
			posStrnds = unbindingComplexSubs(ismembc(unbindingComplexVals, this.complexIndexs_rnaPol2), :);
            if ~isempty(posStrnds)
                %[~, footprint3Prime, footprint5Prime] = this.getDNAFootprint([], this.complexIndexs_rnaPolymerase);
				[~, footprint3Prime, footprint5Prime] = this.getDNAFootprint(i, [], this.complexIndexs_rnaPol2);
                posStrnds( isodd(posStrnds(:, 2)), 1) = posStrnds( isodd(posStrnds(:, 2)), 1) + footprint5Prime;
                posStrnds(iseven(posStrnds(:, 2)), 1) = posStrnds(iseven(posStrnds(:, 2)), 1) + footprint3Prime;
                this.rnaPolymerase.releasePolymerase(posStrnds, proteinIsDegraded);
			end
        end
        function [positionsStrands, sideEffects] = setSiteDamaged(this, i, ...
                damageType, damageSubType, probDamage, maxDamages, ...
                vulnerableMotif, vulnerableMotifType)
				
            import edu.jiangnan.fmme.cell.sim.SimulationStateSideEffect;
            %side effects
            sideEffects = SimulationStateSideEffect.empty(0, 1);
            positionsStrands = zeros(0, 2);
            
            %return if probability of damage is 0
            if probDamage == 0 || maxDamages == 0
                return;
            end
            
            %sample vulnerable sites
            if ischar(vulnerableMotif)
                positionsStrands = this.sampleAccessibleSites(i, probDamage, maxDamages, vulnerableMotif);
            else
                if nnz(this.(vulnerableMotifType){i}) == 0
                    return;
                end
                positionsStrands = find( ...
                    vulnerableMotif == this.(vulnerableMotifType){i} & ...
                    vulnerableMotif == this.damagedSites_nonRedundant{i});
                if isempty(positionsStrands)
                    return;
                end
                maxDamages = min(maxDamages, ...
                    this.randStream.stochasticRound(size(positionsStrands, 1) * probDamage));
                if maxDamages == 0
                    return;
                end
                positionsStrands = this.randStream.randomlySelectNRows(positionsStrands, maxDamages);
            end
            if isempty(positionsStrands)
                return;
            end
            
            %damage selected sites
            this.(damageType)(positionsStrands) = damageSubType;
        end
    end
    %private methods which modify this class' state, and possibly request
    %changes that of other parts of the simulation's state
    methods
        function mergeOwnAdjacentRegions(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            if nnz(this.polymerizedRegions{i}) <= 1
                return;
            end
            
            %% add linking numbers for adjacent double-stranded regions
            
            %convert from SparseMat to coordinates
            [dsPosStrands, lengths]   = find(this.doubleStrandedRegions{i});
            [lkPosStrands, oldLKNums] = find(this.linkingNumbers{i});

            %combine linking numbers for adjacent double-stranded regions
            lkNums = zeros(size(lengths));
            for j = 1:size(lkPosStrands,1)%change from i to j
                idx = find(...
                    lkPosStrands(j, 1) >= dsPosStrands(:, 1)           & ...
                    lkPosStrands(j, 1) <= dsPosStrands(:, 1) + lengths & ...
                    lkPosStrands(j, 2) == dsPosStrands(:, 2));
                lkNums(idx) = lkNums(idx) + oldLKNums(j);
            end
            
            %convert from coordinates to SparseMat
            this.linkingNumbers{i} = CircularSparseMat(...
                dsPosStrands, lkNums, size(this.linkingNumbers{i}), 1);

            %% merge adjacent regions
            this.polymerizedRegions{i} = this.mergeAdjacentRegions(this.polymerizedRegions{i});            
        end
        function polymerizedRegions = mergeAdjacentRegions(~, polymerizedRegions)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            if nnz(polymerizedRegions) <= 1
                return;
            end

            %% merge adjacent regions
            %convert from SparseMat to coordinates
            [positionsStrands, lengths] = find(polymerizedRegions);
            
            %make sure that regions don't overlap
            if any(lengths < 0) || ...
               any(diff(positionsStrands(:,1)) < lengths(1:end-1) & diff(positionsStrands(:,2)) == 0)
                throw(MException('Chromosome:error','polymerizedRegions is corrupt'))
            end
            
            %merge adjacent regions (except over OriC)
            idxs = find(diff(positionsStrands(:,1)) == lengths(1:end-1) & diff(positionsStrands(:,2)) == 0);
            for j = numel(idxs):-1:1%change from i to j
                k(j) = idxs(j); %change from j to k
                lengths(k(j), 1) = lengths(k(j), 1) + lengths(k(j)+1, 1);
            end
            positionsStrands(idxs+1, :) = [];
            lengths(idxs+1, :) = [];

            %% convert from coordinates to SparseMat
            polymerizedRegions = CircularSparseMat(positionsStrands, lengths, size(polymerizedRegions), 1);
        end
        function [rgnPosStrnds, rgnLens] = excludeRegions(this, i, incPosStrnds, incLens, excPosStrnds, excLens)
			L = this.sequenceLen(i);
            %options
            if isscalar(excLens)
                excLens = excLens(ones(size(excPosStrnds,1), 1), 1);
            end
            
            %join included regions
            [incPosStrnds, incLens] = this.joinSplitOverOriCRegions(i, incPosStrnds, incLens);
            
            %exclude excluded regions
            [excPosStrnds, excLens] = this.joinSplitRegions(i, excPosStrnds, excLens);            
            
            excPos = [
                excPosStrnds(:, 1) - L;
                excPosStrnds(:, 1);
                excPosStrnds(:, 1) + L];
            excStrnds = [excPosStrnds(:, 2); excPosStrnds(:, 2); excPosStrnds(:, 2)];
            excLens = [excLens; excLens; excLens];
            
            rgnPos = zeros(0, 1);
            rgnEnds = zeros(0, 1);
            rgnStrnds = zeros(0, 1);
            for j = 1:size(incPosStrnds, 1)
                startCoor = incPosStrnds(j, 1);
                endCoor = startCoor + incLens(j) - 1;
                strnd = incPosStrnds(j, 2);
                
                excIdxs = find(...
                    ((excPos <= startCoor & excPos + excLens-1 >= startCoor) | ...
                    (excPos <= endCoor & excPos + excLens-1 >= endCoor) | ...
                    (excPos >= startCoor & excPos + excLens-1 <= endCoor)) & ...
                    excStrnds == strnd);
                
                if isempty(excIdxs)
                    addtlPos = startCoor;
                    addtlEnds = endCoor;
                elseif excPos(excIdxs(1)) <= startCoor
                    if excPos(excIdxs(end)) + excLens(end) - 1 >= endCoor
                        addtlPos = excPos(excIdxs(1:end-1)) + excLens(excIdxs(1:end-1));
                        addtlEnds = excPos(excIdxs(2:end))-1;
                    else
                        addtlPos = excPos(excIdxs) + excLens(excIdxs);
                        addtlEnds = [excPos(excIdxs(2:end))-1; endCoor];
                    end
                else
                    if excPos(excIdxs(end)) + excLens(end) - 1 >= endCoor
                        addtlPos = [startCoor; excPos(excIdxs(1:end-1)) + excLens(excIdxs(1:end-1))];
                        addtlEnds = excPos(excIdxs)-1;
                    else
                        addtlPos = [startCoor; excPos(excIdxs) + excLens(excIdxs)];
                        addtlEnds = [excPos(excIdxs)-1; endCoor];
                    end
                end
                
                rgnPos = [
                    rgnPos;
                    addtlPos];
                rgnEnds = [
                    rgnEnds;
                    addtlEnds];
                rgnStrnds = [
                    rgnStrnds;
                    strnd(ones(size(addtlPos)), 1)];
            end
            
            idx = find(rgnPos > L);
            rgnPos(idx) = rgnPos(idx) - L;
            rgnEnds(idx) = rgnEnds(idx) - L;
            
            idx = find(rgnPos > rgnEnds);
            rgnPos(idx,:) = [];
            rgnEnds(idx,:) = [];
            rgnStrnds(idx,:) = [];
            
            %join split regions
            [rgnPosStrnds, rgnLens] = this.joinSplitOverOriCRegions(i, [rgnPos rgnStrnds], rgnEnds - rgnPos + 1);
            
            %format output
            rgnPosStrnds(:, 1) = mod(rgnPosStrnds(:, 1) - 1, L) + 1;
            [rgnPosStrnds, order] = edu.jiangnan.fmme.util.SparseMat.sort_subs(rgnPosStrnds, [L this.nCompartments]);
            rgnLens = rgnLens(order, :);
		end
        %finds the intersection of two lists of regions
        function [posStrnds, lens] = intersectRegions(this, i, posStrndsA, lensA, posStrndsB, lensB)
            [posStrndsA lensA] = this.splitOverOriC(i, posStrndsA, lensA);
            posStrndsA(:, 1) = mod(posStrndsA(:, 1) - 1, this.sequenceLen(i)) + 1;
            [posStrndsA, idxs] = edu.jiangnan.fmme.util.SparseMat.sort_subs(posStrndsA, [this.sequenceLen(i) this.nCompartments]);
            lensA = lensA(idxs);
			
            [posStrndsB lensB] = this.splitOverOriC(i, posStrndsB, lensB);  
            posStrndsB(:, 1) = mod(posStrndsB(:, 1) - 1, this.sequenceLen(i)) + 1;
            [posStrndsB, idxs] = edu.jiangnan.fmme.util.SparseMat.sort_subs(posStrndsB, [this.sequenceLen(i) this.nCompartments]);
            lensB = lensB(idxs);

            posStrnds = zeros(0, 2);
            lens = zeros(0, 1);
            for strand = 1:this.nCompartments
                rowsA = posStrndsA(:,2) == strand;
                rowsB = posStrndsB(:,2) == strand;
                posA = posStrndsA(rowsA, 1);
                posB = posStrndsB(rowsB, 1);
                lenA = lensA(rowsA, 1);
                lenB = lensB(rowsB, 1);
                iA = 1;
                iB = 1;
                while iA <= length(posA) && iB <= length(posB)
                    if posA(iA) <= posB(iB)
                        if posA(iA) + lenA(iA) > posB(iB)
                            posStrnds(end+1,:) = [posB(iB) strand];
                            if posA(iA) + lenA(iA) < posB(iB) + lenB(iB)
                                lens(end+1,:) = posA(iA) + lenA(iA) - posB(iB);
                                iA = iA + 1;
                            else
                                lens(end+1,:) = lenB(iB);
                                iB = iB + 1;
                            end
                        else
                            iA = iA + 1;
                        end
                    else
                        if posB(iB) + lenB(iB) > posA(iA)
                            posStrnds(end+1,:) = [posA(iA) strand];
                            if posB(iB) + lenB(iB) < posA(iA) + lenA(iA)
                                lens(end+1,:) = posB(iB) + lenB(iB) - posA(iA);
                                iB = iB + 1;
                            else
                                lens(end+1,:) = lenA(iA);
                                iA = iA + 1;
                            end
                        else
                            iB = iB + 1;
                        end
                    end
                end
            end

            [posStrnds, lens] = this.joinSplitOverOriCRegions(i, posStrnds, lens);
           end
        %split regions into several splitLen sized pieces
        function [posStrnds, lens] = splitRegions(i, ~, rgnPosStrnds, rgnLens, splitLen)
            nPos = sum(floor(abs(rgnLens)/splitLen));
            posStrnds = zeros(nPos, 2);
            lens = repmat(splitLen, nPos, 1);
            
            j = 0;
            for k = 1:size(rgnPosStrnds, 1)% change from i to k
                nPos = floor(abs(rgnLens(k))/splitLen);
                
                posStrnds(j+(1:nPos), 1) = rgnPosStrnds(k, 1) + sign(rgnLens(k))*(0:nPos-1)*splitLen;
                posStrnds(j+(1:nPos), 2) = rgnPosStrnds(k, 2);
                
                j = j + nPos;
            end
        end
        function [posStrnds, lens] = splitOverOriC(this, i, posStrnds, lens)
			idxs = find(posStrnds(:, 1) + lens - 1 > this.sequenceLen(i));
            if isempty(idxs)
                return;
            end
            posStrnds = [posStrnds; ones(numel(idxs), 1) posStrnds(idxs, 2)];
            lens = [lens; posStrnds(idxs, 1) + reshape(lens(idxs), [], 1) - this.sequenceLen(i) - 1];
            lens(idxs) = this.sequenceLen(i) - posStrnds(idxs,1) + 1;
        end
        %join regions which have been split over the ORI
        function [posStrnds, lens] = joinSplitOverOriCRegions(this, i, posStrnds, lens)
		
			pos = posStrnds(:, 1);
            strnds = posStrnds(:, 2);
            ends = pos + lens - 1;
            for j = 1:max(strnds) %change from i to j
                idx1 = find(pos == 1 & strnds == j, 1, 'first');
                idx2 = find(ends == this.sequenceLen(i) & strnds == j, 1, 'first');
                
                if isempty(idx1) || isempty(idx2) || idx1 == idx2
                    continue;
                end
                
                ends(idx2) = ends(idx2) + (ends(idx1) - pos(idx1) + 1);
                pos(idx1, :) = [];
                strnds(idx1, :) = [];
                ends(idx1, :) = [];
            end
            
            posStrnds = [pos strnds];
            lens = ends - pos + 1;
        end
        function [posStrnds, lens] = joinSplitRegions(this, i, posStrnds, lens)
            %sort
            posStrnds(:, 1) = mod(posStrnds(:, 1) - 1, this.sequenceLen(i)) + 1;
            [posStrnds, order] = edu.jiangnan.fmme.util.SparseMat.sort_subs(posStrnds, [this.sequenceLen(i) this.nCompartments]);
            lens = lens(order);
            
            %join
            starts = posStrnds(:, 1);
            ends = starts + lens - 1;
            strnds = posStrnds(:, 2);
            
            tfs = true(size(starts));
            for j = 1:max(strnds) %change from i to j
                idxs = find(strnds == j);
                for k = 1:numel(idxs)-1 %change from j to k
                    if ends(idxs(k))+1 >= starts(idxs(k+1))
                        starts(idxs(k+1)) = starts(idxs(k));
                        ends(idxs(k+1)) = max(ends(idxs(k)), ends(idxs(k+1)));
                        tfs(idxs(k)) = false;
                    end
                end
                if numel(idxs) >= 2
                    idx = idxs(find(tfs(idxs), 1, 'first'));
                    if ends(idxs(end))+1 >= starts(idx)+this.sequenceLen(i)
                        ends(idxs(end)) = this.sequenceLen(i);
                        starts(idx) = 1;
                    end
                end
            end
            
            %format output
            posStrnds = [starts(tfs) strnds(tfs)];
            lens = ends(tfs) - starts(tfs) + 1;
        end
        function [idxs, newIdxs] = excludeOverlappingRegions(this, i, ...
                idxs, newIdxs, positionsStrands, lengths, ...
                footprint, footprint3Prime, footprint5Prime, isPositionsStrandFootprintCentroid, ...
                eitherStrand)
				
            tmpPositionsStrands = positionsStrands([idxs; newIdxs], :);
            tmpLengths = lengths([idxs; newIdxs], 1);
            
            startCoors = tmpPositionsStrands(:, 1) + min(0, tmpLengths + 1);
            if isPositionsStrandFootprintCentroid
                startCoors(mod(tmpPositionsStrands(:, 2), 2) == 1) = ...
                    startCoors(mod(tmpPositionsStrands(:, 2), 2) == 1) - footprint5Prime;
                startCoors(mod(tmpPositionsStrands(:, 2), 2) == 0) = ...
                    startCoors(mod(tmpPositionsStrands(:, 2), 2) == 0) - footprint3Prime;
            end
            endCoors = startCoors + (abs(tmpLengths) - 1) + (footprint - 1);
            if eitherStrand
                strnds = ceil(tmpPositionsStrands(:, 2) / 2);
            else
                strnds = tmpPositionsStrands(:, 2);
            end
            tmpIdxs = [zeros(size(idxs)); (1:numel(newIdxs))'];
            
            tmp = find(startCoors < 0);
            startCoors = [startCoors; startCoors(tmp, :) + this.sequenceLen(i)];
            endCoors   = [endCoors; min(this.sequenceLen(i), endCoors(tmp, :) + this.sequenceLen(i))];
            strnds     = [strnds; strnds(tmp, :)];
            tmpIdxs    = [tmpIdxs; tmpIdxs(tmp, :)];
            
            tmp = find(endCoors > this.sequenceLen(i));
            startCoors = [startCoors; max(1, startCoors(tmp, :) - this.sequenceLen(i))];
            endCoors   = [endCoors; endCoors(tmp, :) - this.sequenceLen(i)];
            strnds     = [strnds; strnds(tmp, :)];
            tmpIdxs    = [tmpIdxs; tmpIdxs(tmp, :)];
            
            tmpTfs = true(size(newIdxs));
            for j = numel(tmpIdxs):-1:2	%change from i to j
                if tmpIdxs(j) <= numel(idxs)
                    continue; 
                end
                tmpTfs(tmpIdxs(j)) = tmpTfs(tmpIdxs(j)) && ~any((...
                    startCoors(1:j-1) <= startCoors(j) & startCoors(j) <= endCoors(1:j-1) | ...
                    startCoors(1:j-1) <= endCoors(j)   & endCoors(j) <= endCoors(1:j-1)) & ...
                    strnds(j) == strnds(1:j-1));
            end
            newIdxs = newIdxs(tmpTfs);
        end
    end
    %setters
    methods        
        %integers [positions x strands] indicating the start positions of
        %polymerized regions of strands and their lengths
        function set.polymerizedRegions(this, value)
            if isequal(this.polymerizedRegions, value)
                return;
            end
            this.polymerizedRegions = value;
            this.validated = this.validated  + 1;
            this.validated_polymerizedRegions = this.validated; %#ok<*MCSUP>
        end
        %integers [positions x strands] indicating the current linking number of
        %each double-stranded region
        function set.linkingNumbers(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.linkingNumbers = value;
            this.validated = this.validated  + 1;
            this.validated_linkingNumbers = this.validated;
        end
        %indices [positions x strands] indicating start positions of protein
        %monomers bound to DNA bases
        function set.monomerBoundSites(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.monomerBoundSites = value;
            this.validated = this.validated  + 1;
            this.validated_proteinBoundSites = this.validated;
        end
        %indices [positions x strands] indicating start positions of
        %macromolecular complexes bound to DNA bases
        function set.complexBoundSites(this, value)
            %NOTE: performance likely better here without checking if new value
            %is different than old
            this.complexBoundSites = value;
            this.validated = this.validated  + 1;
            this.validated_proteinBoundSites = this.validated;
			%end
        end
        %boolean [positions x strands] indicating positions of gap sites
        function set.gapSites(this, value)
            if isequal(this.gapSites, value)
                return;
            end
            this.gapSites = value;
            this.validated = this.validated + 1;
            this.validated_damaged = this.validated;
            this.validated_gapSites = this.validated;
			%end
        end
        %boolean [positions x strands] indicating positions of abasic sites
        function set.abasicSites(this, value)
            if isequal(this.abasicSites, value)
                return;
            end
            this.abasicSites = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_abasicSites = this.validated;
			%end
        end
        %indices [positions x strands] indicating metabolite identity of damaged
        %sugar-phosphates
        function set.damagedSugarPhosphates(this, value)
            if isequal(this.damagedSugarPhosphates, value)
                return;
            end
            this.damagedSugarPhosphates = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_damagedSugarPhosphates = this.validated;
			%end
        end
        %indices [positions x strands] indicating metabolite identity of damaged
        %bases
        function set.damagedBases(this, value)
            if isequal(this.damagedBases, value)
                return;
            end
            this.damagedBases = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_damagedBases = this.validated;
			%end
        end
        %boolean [positions x strands] indicating metabolite identity of
        %intrastrand cross links in DNA
        function set.intrastrandCrossLinks(this, value)
            if isequal(this.intrastrandCrossLinks, value)
                return;
            end
            this.intrastrandCrossLinks = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_intrastrandCrossLinks = this.validated;
			%end
        end
        %boolean [positions x strands] indicating positions of strand breaks in
        %strands of DNA
        function set.strandBreaks(this, value)
            if isequal(this.strandBreaks, value)
                return;
            end
            this.strandBreaks = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_strandBreaks = this.validated;
			%end
        end
        %boolean [positions x strands] indicating positions of holliday
        %junctions
        function set.hollidayJunctions(this, value)
            if isequal(this.hollidayJunctions, value)
                return;
            end
            this.hollidayJunctions = value;
            this.validated = this.validated  + 1;
            this.validated_damaged = this.validated;
            this.validated_hollidayJunctions = this.validated;
			%end
        end
        %boolean indicating whether or not the chromsomes are segregated
        function set.segregated(this, value)
            if isequal(this.segregated, value)
                return;
            end
            this.segregated = value;
            this.validated = this.validated  + 1;
            this.validated_segregated = this.validated;
			end
		%end
    end
    %getters for alternative views of state
    methods
        function invalidate(this)
            %this.validated{i} = uint32(1);
			
			this.validated = ones(16,1);
            
            this.validated_polymerizedRegions     = this.validated;
            this.validated_linkingNumbers         = this.validated;
            this.validated_proteinBoundSites      = this.validated;
            this.validated_damaged                = this.validated;
            this.validated_abasicSites            = this.validated;
            this.validated_gapSites               = this.validated;
            this.validated_damagedSugarPhosphates = this.validated;
            this.validated_damagedBases           = this.validated;
            this.validated_strandBreaks           = this.validated;
            this.validated_intrastrandCrossLinks  = this.validated;
            this.validated_hollidayJunctions      = this.validated;
            this.validated_segregated             = this.validated;
            %{
            this.validated_unpolymerizedRegions          = uint32(0);
            this.validated_singleStrandedRegions         = uint32(0);
            this.validated_doubleStrandedRegions         = uint32(0);
            this.validated_geneCopyNumbers               = uint32(0);
            this.validated_ploidy                        = uint32(0);
            this.validated_polymerizedGenes              = uint32(0);
            this.validated_transcriptionUnitCopyNumbers  = uint32(0);
            this.validated_polymerizedTranscriptionUnits = uint32(0);
            this.validated_geneCopyNumbers_Accessible    = uint32(0);
            this.validated_transcriptionUnitCopyNumbers_Accessible = uint32(0);
            this.validated_accessibleGenes               = uint32(0);
            this.validated_accessibleTranscriptionUnits  = uint32(0);            
            this.validated_linkingNumbers_minFreeEnergy  = uint32(0);
            this.validated_supercoils                    = uint32(0);
            this.validated_superhelicalDensity           = uint32(0);
            this.validated_supercoiled                   = uint32(0);
            this.validated_damagedSites                  = uint32(0);
            this.validated_damagedSites_shifted_incm6AD  = uint32(0);
            this.validated_damagedSites_nonRedundant     = uint32(0);
            this.validated_damagedSites_excm6AD          = uint32(0);
            this.validated_gapSites3                     = uint32(0);
            this.validated_gapSites5                     = uint32(0);
            this.validated_abasicSites3                  = uint32(0);
            this.validated_abasicSites5                  = uint32(0);
            this.validated_damagedSugarPhosphates3       = uint32(0);
            this.validated_damagedSugarPhosphates5       = uint32(0);
            this.validated_damagedBases3                 = uint32(0);
            this.validated_damagedBases5                 = uint32(0);
            this.validated_strandBreaks3                 = uint32(0);
            this.validated_strandBreaks5                 = uint32(0);
            this.validated_intrastrandCrossLinks3        = uint32(0);
            this.validated_intrastrandCrossLinks5        = uint32(0);
            this.validated_hollidayJunctions3            = uint32(0);
            this.validated_hollidayJunctions5            = uint32(0);
            this.validated_singleStrandBreaks            = uint32(0);
            this.validated_doubleStrandBreaks            = uint32(0);
            this.validated_strandBreakClassification     = uint32(0);
            this.validated_munIRMSiteMethylationStatus   = uint32(0);
            this.validated_munIRMSiteRestrictionStatus   = uint32(0);
            this.validated_dryWeight                     = uint32(0);
			%}
			
            this.validated_unpolymerizedRegions          = zeros(16,1);
            this.validated_singleStrandedRegions         = zeros(16,1);
            this.validated_doubleStrandedRegions         = zeros(16,1);
            this.validated_geneCopyNumbers               = zeros(16,1);
            this.validated_ploidy                        = zeros(16,1);
            this.validated_polymerizedGenes              = zeros(16,1);
            this.validated_transcriptionUnitCopyNumbers  = zeros(16,1);
            this.validated_polymerizedTranscriptionUnits = zeros(16,1);
            this.validated_geneCopyNumbers_Accessible    = zeros(16,1);
            this.validated_transcriptionUnitCopyNumbers_Accessible = zeros(16,1);
            this.validated_accessibleGenes               = zeros(16,1);
            this.validated_accessibleTranscriptionUnits  = zeros(16,1);            
            this.validated_linkingNumbers_minFreeEnergy  = zeros(16,1);
            this.validated_supercoils                    = zeros(16,1);
            this.validated_superhelicalDensity           = zeros(16,1);
            this.validated_supercoiled                   = zeros(16,1);
            this.validated_damagedSites                  = zeros(16,1);
            this.validated_damagedSites_shifted_incm6AD  = zeros(16,1);
            this.validated_damagedSites_nonRedundant     = zeros(16,1);
            this.validated_damagedSites_excm6AD          = zeros(16,1);
            this.validated_gapSites3                     = zeros(16,1);
            this.validated_gapSites5                     = zeros(16,1);
            this.validated_abasicSites3                  = zeros(16,1);
            this.validated_abasicSites5                  = zeros(16,1);
            this.validated_damagedSugarPhosphates3       = zeros(16,1);
            this.validated_damagedSugarPhosphates5       = zeros(16,1);
            this.validated_damagedBases3                 = zeros(16,1);
            this.validated_damagedBases5                 = zeros(16,1);
            this.validated_strandBreaks3                 = zeros(16,1);
            this.validated_strandBreaks5                 = zeros(16,1);
            this.validated_intrastrandCrossLinks3        = zeros(16,1);
            this.validated_intrastrandCrossLinks5        = zeros(16,1);
            this.validated_hollidayJunctions3            = zeros(16,1);
            this.validated_hollidayJunctions5            = zeros(16,1);
            this.validated_singleStrandBreaks            = zeros(16,1);
            this.validated_doubleStrandBreaks            = zeros(16,1);
            this.validated_strandBreakClassification     = zeros(16,1);
            this.validated_munIRMSiteMethylationStatus   = zeros(16,1);
            this.validated_munIRMSiteRestrictionStatus   = zeros(16,1);
            this.validated_dryWeight                     = zeros(16,1);			
        end
        
        %integers indicating the start positions of unpolymerized regions (ie.
        %not yet replicated) of strands and their lengths
        function value = get.unpolymerizedRegions(this)
			for i = 1:16
				if this.validated_polymerizedRegions(i) > this.validated_unpolymerizedRegions(i)
					this.unpolymerizedRegions{i} = this.calcUnpolymerizedRegions(i);
					this.validated_unpolymerizedRegions(i) = this.validated(i);
				end
				value = this.unpolymerizedRegions;
			end
		end
			
        function value = calcUnpolymerizedRegions(this, i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
				[polPosStrndsTimes, polLens] = find(this.polymerizedRegions{i});
				polPos = polPosStrndsTimes(:, 1);
				polStrnd = polPosStrndsTimes(:, 2);
				polTimes = polPosStrndsTimes(:, 3:end);
				if isempty(polTimes)
					polTimes = ones(size(polPos));
				end
           
				unpolPos 	= [];
				unpolStrnds = [];
				unpolTimes  = [];
				unpolLens 	= [];
			
				for j = 1:this.nCompartments	%change from i to j
					idxs = find(polStrnd == j);
                
					for k = 1:size(this.polymerizedRegions{i}, 3)	%change from j to k
						idxs2 = idxs(polTimes(idxs) == k);
						if isempty(idxs2)
							unpolPos 	= [unpolPos; 1];
							unpolLens 	= [unpolLens; this.sequenceLen(i)];
							unpolStrnds = [unpolStrnds; j];
							unpolTimes 	= [unpolTimes; k];
							continue;
						end
                    
						unpolPos = [
							unpolPos;
							1;
							polPos(idxs2) + polLens(idxs2)]; %#ok<*AGROW>
						unpolLens = [
							unpolLens;
							polPos(idxs2(1))-1;
							polPos(idxs2(2:end)) - (polPos(idxs2(1:end-1)) + polLens(idxs2(1:end-1)));
							this.sequenceLen(i) - (polPos(idxs2(end)) + polLens(idxs2(end))) + 1;
							]; %#ok<*AGROW>
						unpolStrnds = [
							unpolStrnds;
							j(ones(numel(idxs2) + 1, 1), 1)]; %#ok<*AGROW>
						unpolTimes = [
							unpolTimes;
							k(ones(numel(idxs2)+1, 1), 1)]; %#ok<*AGROW>
					end
				end
            
				idxs = find(unpolLens > 0);
				value = CircularSparseMat([unpolPos(idxs) unpolStrnds(idxs) unpolTimes(idxs)], unpolLens(idxs), size(this.polymerizedRegions{i}), 1);          
        end 
		
        function value = get.singleStrandedRegions(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_singleStrandedRegions(i)
					this.singleStrandedRegions{i} = this.calcSingleStrandedRegions(i);
					this.validated_singleStrandedRegions(i) = this.validated(i);
				end
                
				value = this.singleStrandedRegions;
			end
        end
        function value = calcSingleStrandedRegions(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            [positionsStrandsTimes, lengths] = find(this.polymerizedRegions{i});
			
            strnd = positionsStrandsTimes(:, 2);
            oppStrnd = strnd;
            oppStrnd(mod(strnd, 2) == 0, 1) = strnd(mod(strnd, 2) == 0, 1) - 1;
            oppStrnd(mod(strnd, 2) == 1, 1) = strnd(mod(strnd, 2) == 1, 1) + 1;
            
            starts = mod([
                positionsStrandsTimes(:,1);
                positionsStrandsTimes(:,1) + lengths] ...
                - 1, this.sequenceLen(i)) + 1;
            strandsTimes = [
                strnd positionsStrandsTimes(:, 3:end);
                oppStrnd positionsStrandsTimes(:, 3:end)];
            oppStrandsTimes = [
                oppStrnd positionsStrandsTimes(:, 3:end);
                strnd positionsStrandsTimes(:, 3:end)];
            lengths = [
                lengths;
                repmat(max(lengths), size(lengths))];
			
            [~, ~, ~, extents1] = this.isRegionPolymerized(i, [starts strandsTimes], lengths, true);
			
            [~, ~, ~, extents2] = this.isRegionNotPolymerized(i, [starts oppStrandsTimes], lengths, true);
            
			extents = min(abs(extents1), abs(extents2));
			
            idxs = find(extents > 0);
            
            if size(strandsTimes, 2) == 1
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs(...
					[starts(idxs, :) strandsTimes(idxs, :) extents(idxs)], ...
                    [this.sequenceLen(i) this.nCompartments this.sequenceLen(i)]);
            else
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs(...
					[starts(idxs, :) strandsTimes(idxs, :) extents(idxs)], ...
                    [this.sequenceLen(i) this.nCompartments size(this.polymerizedRegions{i}, 3) this.sequenceLen(i)]);
            end
            
            for j = size(tmp, 1):-1:1
                if tmp(j, 1) == 1
                    idx = find(tmp(:, 1) + tmp(:, 3) - 1 > this.sequenceLen(i) & tmp(:, 2) == tmp(j, 2));                    
                    if ~isempty(idx)
                        tmp(idx, 3) = max(tmp(idx, 3), (this.sequenceLen(i) - tmp(idx,1) + 1) + tmp(j, 3));
                        tmp(j, :) = [];
                    end
                end
            end
            
            posStrnds = tmp(:, 1:end-1);
            lens = tmp(:, end);
            [posStrnds, lens] = this.splitOverOriC(i, posStrnds, lens);
            value = this.mergeAdjacentRegions(CircularSparseMat(posStrnds, lens, size(this.polymerizedRegions{i}), 1));
        end
                
        function value = get.doubleStrandedRegions(this)
			for i = 1:16
            if  this.validated_polymerizedRegions(i) > this.validated_doubleStrandedRegions(i)
                this.doubleStrandedRegions{i} = this.calcDoubleStrandedRegions(i);
                this.validated_doubleStrandedRegions(i) = this.validated(i);
            end
                
            value = this.doubleStrandedRegions;
			end
        end
		
        function value = calcDoubleStrandedRegions(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            [positionsStrandsTimes, lengths] = find(this.polymerizedRegions{i});
            
            if size(positionsStrandsTimes, 2) == 2
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs([
                    positionsStrandsTimes(:, 1) ceil(positionsStrandsTimes(:, 2)/2) lengths;
                    ], [this.sequenceLen(i) this.nCompartments this.sequenceLen(i)]); 
				
                [~, ~, ~, extents1] = this.isRegionPolymerized(i, [tmp(:, 1) 2*tmp(:,2)-1], tmp(:, end), true);
                [~, ~, ~, extents2] = this.isRegionPolymerized(i, [tmp(:, 1) 2*tmp(:,2)  ], tmp(:, end), true);
                extents = min(abs(extents1), abs(extents2));
                idxs = find(extents > 0);
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs([
                    tmp(idxs, 1) 2*tmp(idxs, 2)-1 extents(idxs)
                    tmp(idxs, 1) 2*tmp(idxs, 2)   extents(idxs)
                    ], [this.sequenceLen(i) this.nCompartments this.sequenceLen(i)]);
            else
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs([
                    positionsStrandsTimes(:, 1) ceil(positionsStrandsTimes(:, 2)/2) positionsStrandsTimes(:, 3:end) lengths;
                    ], [this.sequenceLen(i) this.nCompartments size(this.polymerizedRegions, 3) this.sequenceLen(i)]);
                
                [~, ~, ~, extents1] = this.isRegionPolymerized(i, [tmp(:, 1) 2*tmp(:,2)-1 tmp(:,3:end-1)], tmp(:, end), true);
                [~, ~, ~, extents2] = this.isRegionPolymerized(i, [tmp(:, 1) 2*tmp(:,2)   tmp(:,3:end-1)], tmp(:, end), true);
                extents = min(abs(extents1), abs(extents2));
                idxs = find(extents > 0);
                
                tmp = edu.jiangnan.fmme.util.SparseMat.unique_subs([
                    tmp(idxs, 1) 2*tmp(idxs, 2)-1 tmp(idxs, 3:end-1) extents(idxs)
                    tmp(idxs, 1) 2*tmp(idxs, 2)   tmp(idxs, 3:end-1) extents(idxs)
                    ], [this.sequenceLen(i) this.nCompartments size(this.polymerizedRegions{i}, 3) this.sequenceLen(i)]);
            end
            
            for j = size(tmp, 1):-1:1
                if tmp(j, 1) == 1
                    idx = find(tmp(:, 1) + tmp(:, 3) - 1 > this.sequenceLen(i) & tmp(:, 2) == tmp(j, 2));
                    if ~isempty(idx)
                        tmp(idx, 3) = max(tmp(idx, 3), (this.sequenceLen(i) - tmp(idx,1) + 1) + tmp(j, 3));
                        tmp(j, :) = [];
                    end
                end
            end
            
            posStrnds = tmp(:, 1:end-1);
            lens = tmp(:, end);
            [posStrnds, lens] = this.splitOverOriC(i, posStrnds, lens);
            value = this.mergeAdjacentRegions(CircularSparseMat(posStrnds, lens, size(this.polymerizedRegions{i}), 1));
        end
        
        function value = get.geneCopyNumbers(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_geneCopyNumbers(i)
					this.geneCopyNumbers{i} = this.calcGeneCopyNumbers(i);
					this.validated_geneCopyNumbers(i) = this.validated(i);
				end
                
				value = this.geneCopyNumbers;
			end
        end
		
        %number of copies of each gene that have been polymerized (Nx1)
        function value = calcGeneCopyNumbers(this,i)
            value = sum(this.polymerizedGenes{i}, 2);
        end
		
        function value = get.ploidy(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_ploidy(i)
					this.ploidy{i} = this.calcPloidy(i);
					this.validated_ploidy(i) = this.validated(i);
				end
            
				value = this.ploidy;
			end
        end
		
        function value = calcPloidy(this,i)
		
            value = collapse(this.polymerizedRegions{i})/(2*this.sequenceLen(i));
			
        end
		
        function value = get.polymerizedGenes(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_polymerizedGenes(i)
					this.polymerizedGenes{i} = this.calcPolymerizedGenes(i);
					this.validated_polymerizedGenes(i) = this.validated(i);
				end
                
				value = this.polymerizedGenes;
			end
		end
		
        %whether each copy of each gene has been polymerized (Nx2)
        function value = calcPolymerizedGenes(this,i)
			value = this.isRegionPolymerized(i, ...
                [this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i};
                 this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i}+2],...
                [this.transcriptionUnitLengths{i}; this.transcriptionUnitLengths{i}], false);
            
			value = reshape(value, [], 2);	
        end
        
        
        function value = get.transcriptionUnitCopyNumbers(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_transcriptionUnitCopyNumbers(i)
					this.transcriptionUnitCopyNumbers{i} = this.getTranscriptionUnitCopyNumbers(i);
					this.validated_transcriptionUnitCopyNumbers(i) = this.validated(i);
				end
                
				value = this.transcriptionUnitCopyNumbers;
			end
		end
		
        %number of copies of each transcription unit that have been polymerized (Nx1)
        function value = getTranscriptionUnitCopyNumbers(this,i)
		
            value = sum(this.polymerizedTranscriptionUnits{i}, 2);
        end

        
        function value = get.polymerizedTranscriptionUnits(this)
			for i = 1:16
				if  this.validated_polymerizedRegions(i) > this.validated_polymerizedTranscriptionUnits(i)
					this.polymerizedTranscriptionUnits{i} = this.calcPolymerizedTranscriptionUnits(i);
					this.validated_polymerizedTranscriptionUnits(i) = this.validated(i);
				end
                
				value = this.polymerizedTranscriptionUnits;
			end
        end
        %whether each copy of each transcription unit has been polymerized (Nx2)
        function value = calcPolymerizedTranscriptionUnits(this,i)
            value = this.isRegionPolymerized(i, ...
                [this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i};
                 this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i}+2],...
                [this.transcriptionUnitLengths{i}; this.transcriptionUnitLengths{i}], false);
            
			value = reshape(value, [], 2);
        end
        
        function value = get.geneCopyNumbers_Accessible(this)
			for i = 1:16
				val{i} = max([
                this.validated_polymerizedRegions(i); 
                this.validated_damaged(i);
                this.validated_proteinBoundSites(i)]);
            if val{i} > this.validated_geneCopyNumbers_Accessible(i)
                this.geneCopyNumbers_Accessible{i} = this.calcCopyNumbers_Accessible(i);
                this.validated_geneCopyNumbers_Accessible(i) = this.validated(i);
            end
                
            value = this.geneCopyNumbers_Accessible;
			end
        end

        %number of copies of each gene that are accessible
        function value = calcCopyNumbers_Accessible(this, i)
            value = sum(this.accessibleGenes{i}, 2);
        end
        
        
        function value = get.transcriptionUnitCopyNumbers_Accessible(this)
			for i = 1:16
            val{i} = max([
                this.validated_polymerizedRegions(i); 
                this.validated_damaged(i);
                this.validated_proteinBoundSites(i)]);
            if val{i} > this.validated_transcriptionUnitCopyNumbers_Accessible(i)
                this.transcriptionUnitCopyNumbers_Accessible{i} = this.calcTranscriptionUnitCopyNumbers_Accessible(i);
                this.validated_transcriptionUnitCopyNumbers_Accessible(i) = this.validated(i);
            end
                
            value = this.transcriptionUnitCopyNumbers_Accessible;
			end
        end
                
        %number of copies of each transcription unit that are accessible
        function value = calcTranscriptionUnitCopyNumbers_Accessible(this,i)
            value = sum(this.accessibleTranscriptionUnits{i}, 2);
        end
        
        
        function value = get.accessibleGenes(this)
			for i = 1:16
				val{i} = max([
					this.validated_polymerizedRegions(i); 
					this.validated_damaged(i);
					this.validated_proteinBoundSites(i)]);
            if 	val{i} > this.validated_accessibleGenes(i)
                this.accessibleGenes{i} = this.calcAccessibleGenes(i);
                this.validated_accessibleGenes(i) = this.validated(i);
            end
                
            value = this.accessibleGenes;
			end
        end
        
        %boolean indicator of undamaged, unoccupied genes
        %true  ==> gene is accessible
        %false ==> gene is inaccessible
        function value = calcAccessibleGenes(this,i)
            value = reshape(this.isRegionAccessible(i, [...
				this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i};
                this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i}+2],...
                [this.transcriptionUnitLengths{i}; this.transcriptionUnitLengths{i}], [], [], ...
				true, [], false, true), [], this.nCompartments/2);
        end
        
        
        function value = get.accessibleTranscriptionUnits(this)
			for i = 1: 16
				val{i} = max([
					this.validated_polymerizedRegions(i); 
					this.validated_damaged(i);
					this.validated_proteinBoundSites(i)]);
				if val{i} > this.validated_accessibleTranscriptionUnits(i)
					this.accessibleTranscriptionUnits{i} = this.calcAccessibleTranscriptionUnits(i);
					this.validated_accessibleTranscriptionUnits(i) = this.validated(i);
				end
                
				value = this.accessibleTranscriptionUnits;
			end
        end
        
        %boolean indicator of undamaged, unoccupied transcription units
        %true  ==> transcription unit is accessible
        %false ==> transcription unit is inaccessible
        function value = calcAccessibleTranscriptionUnits(this,i)
            value = reshape(this.isRegionAccessible(i, [...
                this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i};
                this.transcriptionUnitStartCoordinates{i} this.transcriptionUnitStrands{i}+2], ...
                [this.transcriptionUnitLengths{i}; this.transcriptionUnitLengths{i}], [], [], true, [], false, true), ...
                [], this.nCompartments/2);
        end
        
        function value = get.strandBreakClassification(this)
			for i = 1: 16
				if max(this.validated_polymerizedRegions(i), this.validated_damaged(i)) > this.validated_strandBreakClassification(i)
					this.strandBreakClassification{i} = this.calcStrandBreakClassification(i);
					this.validated_strandBreakClassification(i) = this.validated(i);
				end
                
				value = this.strandBreakClassification;
			end
        end
        %numbers of each class of strand break (SSB, SSB+, 2SSB, DSB, DSB+,
        %DSB++) in DNA
        function value = calcStrandBreakClassification(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            import edu.jiangnan.fmme.util.countUnique;
            %parameters
            segmentLength = this.strandBreakClassification_segmentLength;
            dsbSep = this.strandBreakClassification_doubleStrandBreakSeparation;
            genomeLength = this.sequenceLen(i);
            numStrands = this.nCompartments;
            numTime = size(this.strandBreaks, 3);
            numSegments = ceil(genomeLength / segmentLength);
            
            %damaged Sites
            [polymerizedPositionsStrands, polymerizedLengths] = find(this.polymerizedRegions{i});
            idxs = find(polymerizedPositionsStrands(:,1) <= segmentLength * numSegments - genomeLength);
            polymerizedPositionsStrands = [polymerizedPositionsStrands; polymerizedPositionsStrands(idxs, 1) + this.sequenceLen(i) polymerizedPositionsStrands(idxs, 2:end)];
            polymerizedLengths = [polymerizedLengths; polymerizedLengths(idxs, :)];
            idxs = find(polymerizedPositionsStrands(:,1) + polymerizedLengths - 1 > segmentLength * numSegments);
            polymerizedLengths(idxs,:) = segmentLength * numSegments - polymerizedPositionsStrands(idxs, 1) + 1;
            [polymerizedPositionsStrands, polymerizedLengths] = find(this.mergeAdjacentRegions(...
                CircularSparseMat(polymerizedPositionsStrands, polymerizedLengths, [segmentLength * numSegments numStrands numTime], 1)));
            
            if numTime == 1
                polymerizedPositionsStrands = [polymerizedPositionsStrands ones(size(polymerizedPositionsStrands,1), 1)];
            end
            
            polymerizedStrands = [];
            for j = 1:size(polymerizedPositionsStrands, 1)%change from i to j
                regions = ceil((polymerizedPositionsStrands(j,1)-1) / segmentLength)+1 : floor((polymerizedPositionsStrands(j,1) + polymerizedLengths(j,1) -1) / segmentLength);
                polymerizedStrands = [polymerizedStrands;
                    regions' ...
                    repmat(ceil(polymerizedPositionsStrands(j,2)/2), numel(regions), 1) ...
                    repmat(polymerizedPositionsStrands(j,3), numel(regions), 1)];
            end
            [polymerizedRegions, ~, idxs]= unique(polymerizedStrands, 'rows'); %#ok<PROP>
            [idxs, counts] = countUnique(idxs);
            polymerizeTimes = polymerizedRegions(idxs(counts == 2), 3:end); %#ok<PROP>
            
            numPolymerizedSegments = zeros(1, 1, numTime);
            [idxs, counts] = countUnique(polymerizeTimes);
            numPolymerizedSegments(idxs) = counts;
            
            %initialize classification
            dmgSites = [this.strandBreaks{i}; this.strandBreaks{i}(1:segmentLength * numSegments - genomeLength, :, :)];
            
            subs = find(permute(reshape(dmgSites, [], numSegments, numStrands, numTime), [3 1 2 4])); %[strands X positions X segments X time]
            if numTime == 1; subs = [subs ones(size(subs,1), 1)]; end
            segmentTimeInds = sub2ind([numSegments numTime], subs(:, 3), subs(:, 4));
            
            value = zeros(this.strandBreakClassification_index_DSB__, 1, numTime);
            
            %classify segments
            while ~isempty(subs)
                %time
                time = subs(1, 4);
                
                %find extent of segment
                endIdx = find(...
                    segmentTimeInds(1) ~= segmentTimeInds | ...
                    ceil(subs(1,1)/2) ~= ceil(subs(:,1)/2), ...
                    1, 'first') - 1;
                if isempty(endIdx)
                    endIdx = size(subs, 1);
                end
                
                if ~this.isRegionDoubleStranded(i, [(subs(:,3)-1)*segmentLength(i)+1 subs(:,1) subs(:,4)], segmentLength, false);
                    continue;
                end
                
                %classify damaged segment
                if endIdx == 1
                    classification = this.strandBreakClassification_index_SSB;
                else
                    strands   = subs(1:endIdx, 1);
                    positions = subs(1:endIdx, 2);
                    if ~isempty(strands) && all(strands == strands(1))
                        classification = this.strandBreakClassification_index_SSB_;
                    else
                        dsbIdx = find(diff(positions) < dsbSep & diff(strands));
                        if isempty(dsbIdx)
                            classification = this.strandBreakClassification_index_2SSB;
                        else
                            if numel(dsbIdx)>1 && diff(positions(dsbIdx([1 end]))) >= dsbSep
                                classification = this.strandBreakClassification_index_DSB__;
                            elseif numel(dsbIdx)>1 || ...
                                    (dsbIdx>1 && diff(positions([dsbIdx - 1 dsbIdx])) < dsbSep) || ...
                                    (dsbIdx<numel(positions) - 1 && diff(positions([dsbIdx + 1 dsbIdx + 2])) < dsbSep)
                                classification = this.strandBreakClassification_index_DSB_;
                            else
                                classification = this.strandBreakClassification_index_DSB;
                            end
                        end
                    end
                end
                
                %update counts of classified segments
                value(classification, 1, time) = value(classification, 1, time) + 1;
                
                %shrink subs
                subs(1:endIdx, :) = [];
                segmentTimeInds(1:endIdx, :) = [];
            end
            
            %compute numbers of segments without damage
            value(this.strandBreakClassification_index_NB, 1, :) =  numPolymerizedSegments - sum(value,1);
        end
		
        function value = get.linkingNumbers_minFreeEnergy(this)
			for i = 1: 16
            if this.validated_polymerizedRegions(i) > this.validated_linkingNumbers_minFreeEnergy(i)
                this.linkingNumbers_minFreeEnergy{i} = this.calcLinkingNumbers_minFreeEnergy(i);
                this.validated_linkingNumbers_minFreeEnergy(i) = this.validated(i);
            end
            
            value = this.linkingNumbers_minFreeEnergy;
        end
        end
		
        function value = calcLinkingNumbers_minFreeEnergy(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            [posStrands, lens] = find(this.doubleStrandedRegions{i});
            
            value = CircularSparseMat(posStrands, lens / this.relaxedBasesPerTurn(i), size(this.linkingNumbers{i}), 1);
        end
		
        function value = get.supercoils(this)
			for i = 1: 16
				if  any([this.validated_polymerizedRegions(i); this.validated_linkingNumbers(i)] > this.validated_supercoils(i))
					this.supercoils{i} = this.calcSupercoils(i);
					this.validated_supercoils(i) = this.validated(i);
				end
            
				value = this.supercoils;
			end
        end
		
        function value = calcSupercoils(this,i)
		
            value = this.linkingNumbers{i} - this.linkingNumbers_minFreeEnergy{i};
        end
		
        function value = get.superhelicalDensity(this)
			for i = 1: 16
				if  any([this.validated_polymerizedRegions(i); this.validated_linkingNumbers(i)] > this.validated_superhelicalDensity(i))
					this.superhelicalDensity{i} = this.calcSuperhelicalDensity(i);
					this.validated_superhelicalDensity(i) = this.validated(i);
				end
            
            value = this.superhelicalDensity;
			end
        end
		
        function value = calcSuperhelicalDensity(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            lkNums_min = this.linkingNumbers_minFreeEnergy{i};
            [posStrnds, deltas] = find(this.linkingNumbers{i} - lkNums_min);
			
            value = CircularSparseMat(posStrnds, deltas ./ lkNums_min(posStrnds), size(this.linkingNumbers{i}), 1);
        end
		
        %check if superhelical density within tolerance of equilbrium value
        function value = get.supercoiled(this)
			for i = 1: 16
				if any([this.validated_polymerizedRegions(i); this.validated_linkingNumbers(i)] > this.validated_supercoiled(i))
					this.supercoiled{i} = this.calcSupercoiled(i);
					this.validated_supercoiled(i) = this.validated(i);
				end
            
				value = this.supercoiled;
			end
        end
		
        function value = calcSupercoiled(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            siz = [this.sequenceLen(i) this.nCompartments];
            
            [posStrnds, lens] = find(this.doubleStrandedRegions{i});
            lks = this.linkingNumbers{i}(posStrnds);
            
            lk0s = lens / this.relaxedBasesPerTurn(i);
            sigmas = (lks - lk0s) ./ lk0s;
            tfs = abs(sigmas - this.equilibriumSuperhelicalDensity) < this.supercoiledSuperhelicalDensityTolerance;
            
            value = CircularSparseMat(posStrnds, tfs, siz, 1);
        end
		
        function value = get.damagedSites(this)
			for i = 1:16
				if  this.validated_damaged(i) > this.validated_damagedSites(i)
					this.damagedSites{i} = this.calcDamagedSites(i);
					this.validated_damagedSites(i) = this.validated(i);
				end
            
				value = this.damagedSites;
			end
        end
		
        function value = calcDamagedSites(this,i)
            value = this.getDamagedSites(i, true, true, true, true, false, true, false);
        end
        
        function value = get.damagedSites_shifted_incm6AD(this)
			for i = 1:16
				if  this.validated_damaged(i) > this.validated_damagedSites_shifted_incm6AD(i)
					this.damagedSites_shifted_incm6AD{i} = this.calcDamagedSites_shifted_incm6AD(i);
					this.validated_damagedSites_shifted_incm6AD(i) = this.validated(i);
            end
            
				value = this.damagedSites_shifted_incm6AD;
			end
        end
        
        function value = calcDamagedSites_shifted_incm6AD(this,i)
            value = this.getDamagedSites(i, true, true, true, true, false, true, true);
        end
        
        function value = get.damagedSites_nonRedundant(this)
			for i = 1:16
				if  this.validated_damaged(i) > this.validated_damagedSites_nonRedundant(i)
					this.damagedSites_nonRedundant{i} = this.calcDamagedSites_nonRedundant(i);
					this.validated_damagedSites_nonRedundant(i) = this.validated(i);
				end
            
				value = this.damagedSites_nonRedundant;
			end
        end
        
        function value = calcDamagedSites_nonRedundant(this,i)
            value = this.getDamagedSites(i, true, true, false, false, false, false, true);
        end
        
        function value = get.damagedSites_excm6AD(this)
			for i = 1:16
				if  this.validated_damaged(i) > this.validated_damagedSites_excm6AD(i)
					this.damagedSites_excm6AD{i} = this.calcDamagedSites_excm6AD(i);
					this.validated_damagedSites_excm6AD(i) = this.validated(i);
				end
            
				value = this.damagedSites_excm6AD;
			end
        end
        
        function value = calcDamagedSites_excm6AD(this,i)
            value = this.getDamagedSites(i, true, true, false, false, false, false, false);
        end
        
        %boolean (genome length x 2) indicating positions of gap sites 3'
        %to bases
        function value = get.gapSites3(this)
			for i = 1:16
				if  this.validated_gapSites(i) > this.validated_gapSites3(i)
					this.gapSites3{i} = this.calcGapSites3(i);
					this.validated_gapSites3(i) = this.validated(i);
				end
            
				value = this.gapSites3;
			end
        end
        
        function value = calcGapSites3(this,i)
            value = this.shiftCircularSparseMatBase5Prime(this.gapSites{i});
        end
        
        %boolean (genome length x 2) indicating positions of gap sites 5'
        %to bases
        function value = get.gapSites5(this)
			for i = 1:16
				if  this.validated_gapSites(i) > this.validated_gapSites5(i)
					this.gapSites5{i} = this.calcGapSites5(i);
					this.validated_gapSites5(i) = this.validated(i);
            end
            
				value = this.gapSites5;
			end
        end
        
        function value = calcGapSites5(this,i)
            value = this.shiftCircularSparseMatBase3Prime(this.gapSites{i});
        end
        
        %boolean (genome length x 2) indicating positions of abasic sites 3'
        %to bases
        function value = get.abasicSites3(this)
			for i = 1:16
				if  this.validated_abasicSites(i) > this.validated_abasicSites3(i)
					this.abasicSites3{i} = this.calcAbasicSites3(i);
					this.validated_abasicSites3(i) = this.validated(i);
            end
            
				value = this.abasicSites3;
			end
        end
        
        function value = calcAbasicSites3(this,i)
            value = this.shiftCircularSparseMatBase5Prime(this.abasicSites{i});
        end
        
        %boolean (genome length x 2) indicating positions of abasic sites 5'
        %to bases
        function value = get.abasicSites5(this)
			for i = 1:16
				if  this.validated_abasicSites(i) > this.validated_abasicSites5(i)
					this.abasicSites5{i} = this.calcAbasicSites5(i);
					this.validated_abasicSites5(i) = this.validated(i);
				end
            
				value = this.abasicSites5;
			end
        end
        
        function value = calcAbasicSites5(this,i)
            value = this.shiftCircularSparseMatBase3Prime(this.abasicSites{i});
        end
        
        %integer (genome length x 2) indicating indices of damaged sugar
        %phosphates 3' to bases
        function value = get.damagedSugarPhosphates3(this)
			for i = 1:16
				if this.validated_damagedSugarPhosphates(i) > this.validated_damagedSugarPhosphates3(i)
					this.damagedSugarPhosphates3{i} = this.calcDamagedSugarPhosphates3(i);
					this.validated_damagedSugarPhosphates3(i) = this.validated(i);
				end
            
				value = this.damagedSugarPhosphates3;
			end
        end
        
        function value = calcDamagedSugarPhosphates3(this,i)
            value = this.shiftCircularSparseMatBase5Prime(this.damagedSugarPhosphates{i});
        end
        
        %integer (genome length x 2) indicating indices of damaged sugar
        %phosphates 5' to bases
        function value = get.damagedSugarPhosphates5(this)
			for i = 1:16
				if  this.validated_damagedSugarPhosphates(i) > this.validated_damagedSugarPhosphates5(i)
					this.damagedSugarPhosphates5{i} = this.calcDamagedSugarPhosphates5(i);
					this.validated_damagedSugarPhosphates5(i) = this.validated(i);
				end
            
				value = this.damagedSugarPhosphates5;
			end
        end
        
        function value = calcDamagedSugarPhosphates5(this,i)
            value = this.shiftCircularSparseMatBase3Prime(this.damagedSugarPhosphates{i});
        end
        
        %integer (genome length x 2) indicating indices of damaged bases 3' to bases
        function value = get.damagedBases3(this)
			for i = 1:16
				if this.validated_damagedBases(i) > this.validated_damagedBases3(i)
                this.damagedBases3{i} = this.calcDamagedBases3(i);
                this.validated_damagedBases3(i) = this.validated(i);
				end
            
				value = this.damagedBases3;
			end
        end
        
        function value = calcDamagedBases3(this,i)
            value = this.shiftCircularSparseMatBase5Prime(this.damagedBases{i});
        end
        
        %integer (genome length x 2) indicating indices of damaged bases 5' to bases
        function value = get.damagedBases5(this)
			for i = 1:16
				if this.validated_damagedBases(i) > this.validated_damagedBases5(i)
                this.damagedBases5{i} = this.calcDamagedBases5(i);
                this.validated_damagedBases5(i) = this.validated(i);
				end
            
				value = this.damagedBases5;
			end
        end
        
        function value = calcDamagedBases5(this,i)
            value = this.shiftCircularSparseMatBase3Prime(this.damagedBases{i});
        end
        
        %boolean (genome length x 2) indicating positions of intrastrand cross
        %links 3' to bases
        function value = get.intrastrandCrossLinks3(this)
			for i = 1:16
				if  this.validated_intrastrandCrossLinks(i) > this.validated_intrastrandCrossLinks3(i)
					this.intrastrandCrossLinks3{i} = this.calcIntrastrandCrossLinks3(i);
					this.validated_intrastrandCrossLinks3(i) = this.validated(i);
				end
            
				value = this.intrastrandCrossLinks3;
			end
        end
        
        function value = calcIntrastrandCrossLinks3(this,i)
            value = this.shiftCircularSparseMatBase5Prime(this.intrastrandCrossLinks{i});
        end
        
        %boolean (genome length x 2) indicating positions of intrastrand cross
        %links 5' to bases
        function value = get.intrastrandCrossLinks5(this)
			for i = 1:16
				if  this.validated_intrastrandCrossLinks(i) > this.validated_intrastrandCrossLinks5(i)
					this.intrastrandCrossLinks5{i} = this.calcIntrastrandCrossLinks5(i);
					this.validated_intrastrandCrossLinks5(i) = this.validated(i);
            end
            
				value = this.intrastrandCrossLinks5;
			end
        end
        
        function value = calcIntrastrandCrossLinks5(this,i)
            value = this.shiftCircularSparseMatBase3Prime(this.intrastrandCrossLinks{i});
        end
        
        %boolean (genome length x 2) indicating positions of strand breaks
        %3' to bases
        function value = get.strandBreaks3(this)
			for i = 1:16
				if  this.validated_strandBreaks(i) > this.validated_strandBreaks3(i)
					this.strandBreaks3{i} = this.calcStrandBreaks3(i);
					this.validated_strandBreaks3(i) = this.validated(i);
				end
            
				value = this.strandBreaks3;
			end
        end
        
        function value = calcStrandBreaks3(this,i)
            value = this.unshiftCircularSparseMatBond3Prime(this.strandBreaks{i});
        end
        
        %boolean (genome length x 2) indicating positions of strand breaks
        %5' to bases
        function value = get.strandBreaks5(this)
			for i = 1:16
				if this.validated_strandBreaks(i) > this.validated_strandBreaks5(i)
					this.strandBreaks5{i} = this.calcStrandBreaks5(i);
					this.validated_strandBreaks5(i) = this.validated(i);
				end
            
				value = this.strandBreaks5;
			end
        end
        
        function value = calcStrandBreaks5(this,i)
            value = this.unshiftCircularSparseMatBond5Prime(this.strandBreaks{i});
        end
        
        %boolean (genome length x 2) indicating positions of holliday junctions
        %3' to bases
        function value = get.hollidayJunctions3(this)
			for i = 1:16
				if  this.validated_hollidayJunctions(i) > this.validated_hollidayJunctions3(i)
					this.hollidayJunctions3{i} = this.calcHollidayJunctions3(i);
					this.validated_hollidayJunctions3(i) = this.validated(i);
				end
            
				value = this.hollidayJunctions3;
			end
        end
        
        function value = calcHollidayJunctions3(this,i)
            value = this.shiftCircularSparseMatBond5Prime(this.hollidayJunctions{i});
        end
        
        %boolean (genome length x 2) indicating positions of holliday junctions
        %5' to bases
        function value = get.hollidayJunctions5(this)
			for i = 1:16
            if this.validated_hollidayJunctions(i) > this.validated_hollidayJunctions5(i)
                this.hollidayJunctions5{i} = this.calcHollidayJunctions5(i);
                this.validated_hollidayJunctions5(i) = this.validated(i);
            end
            
            value = this.hollidayJunctions5;
			end
        end
        
        function value = calcHollidayJunctions5(this,i)
            value = this.shiftCircularSparseMatBond3Prime(this.hollidayJunctions{i});
        end
        
        %boolean (genome length x 2) indicating positions of single
        %strand breaks -- strand breaks excluding
        %- double strand breaks that are part of double strand breaks
        %- strand breaks adjacent to gap sites
        function value = get.singleStrandBreaks(this)
			for i = 1:16
            if this.validated_damaged(i) > this.validated_singleStrandBreaks(i)
                this.singleStrandBreaks{i} = this.calcSingleStrandBreaks(i);
                this.validated_singleStrandBreaks(i) = this.validated(i);
            end
            
            value = this.singleStrandBreaks;
			end
        end
        
        function value = calcSingleStrandBreaks(this,i)
            value = this.strandBreaks{i};
            
            %exclude strand breaks that are part of double strand breaks
            value(find(this.doubleStrandBreaks{i})) = 0; %#ok<FNDSB>
            
            %exclude strand breaks adjacent to other damage (except holliday
            %junctions)
            otherDamages = this.getDamagedSites(i, true, false, true, true, false, false, false);
            value(find( ...
                this.shiftCircularSparseMatBond3Prime(otherDamages) | ...
                this.shiftCircularSparseMatBond5Prime(otherDamages) ...
                )) = 0; %#ok<FNDSB>
            
            %cast to logical sparse mat
            value = valueCast(value, 'logical');
        end
        
        function value = get.doubleStrandBreaks(this)
			for i = 1:16
            if this.validated_strandBreaks(i) > this.validated_doubleStrandBreaks(i)
                this.doubleStrandBreaks{i} = this.calcDoubleStrandBreaks(i);
                this.validated_doubleStrandBreaks(i) = this.validated(i);
            end
            
            value = this.doubleStrandBreaks;
			end
        end
        
        function value = calcDoubleStrandBreaks(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if this.doubleStrandBreakSeparation ~= 1
                throw(MException('DNARepair:error','Simulation only valid for doubleStrandBreakSeparation equal 1'));
            end
            
            value = this.strandBreaks{i}(:, 1:2:end, :) & this.strandBreaks{i}(:, 2:2:end, :);
            value = [value(:,1) value(:,1) value(:,2) value(:,2)];
        end
                
        function value = get.restrictableMunIRMSites(this)
			for i = 1:16
				if ...
                    this.validated_damagedBases(i) > this.validated_munIRMSiteRestrictionStatus(i) || ...
                    this.validated_strandBreaks(i) > this.validated_munIRMSiteRestrictionStatus(i)
                this.restrictableMunIRMSites{i} = this.calcRestrictableMunIRMSites(i);
                this.validated_munIRMSiteRestrictionStatus(i) = this.validated(i);
				end
            
				value = this.restrictableMunIRMSites;
			end
        end
        
        function value = get.hemiunmethylatedMunIRMSites(this)
			for i = 1:16
				if  this.validated_damagedBases(i) > this.validated_munIRMSiteMethylationStatus(i)
					this.hemiunmethylatedMunIRMSites{i} = this.calcHemiunmethylatedMunIRMSites(i);
					this.validated_munIRMSiteMethylationStatus(i) = this.validated(i);
				end
            
				value = this.hemiunmethylatedMunIRMSites;
			end
        end
        
        function value = calcRestrictableMunIRMSites(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            dr = this.dnaRepair;
            nSites = size(dr.RM_MunI_RecognitionSites{i}, 1);
            methylationPosStrnds = [
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(2)) 4*ones(nSites, 1)];
            restrictionPosStrnds = [
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(2)) 4*ones(nSites, 1)];
            
            isPositionMethylated = reshape(this.damagedBases{i}(methylationPosStrnds), [], 2) == this.metabolite.m6ADIndexs;
			
			isPositionDamaged = [
                cat(3, ...
                    reshape(this.damagedSites_shifted_incm6AD{i}([dr.RM_MunI_RecognitionSites{i}(:)   ones(6*nSites, 1)]), [], 6), ...
                    reshape(this.damagedSites_shifted_incm6AD{i}([dr.RM_MunI_RecognitionSites{i}(:) 2*ones(6*nSites, 1)]), [], 6))
                cat(3, ...
                    reshape(this.damagedSites_shifted_incm6AD{i}([dr.RM_MunI_RecognitionSites{i}(:) 3*ones(6*nSites, 1)]), [], 6), ...
                    reshape(this.damagedSites_shifted_incm6AD{i}([dr.RM_MunI_RecognitionSites{i}(:) 4*ones(6*nSites, 1)]), [], 6))];
            isPositionStrandBreaks = [
                cat(3, ...
                    this.strandBreaks{i}([reshape(dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(1)), [], 1),   ones(nSites, 1)]),...
                    this.strandBreaks{i}([reshape(dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(2)), [], 1), 2*ones(nSites, 1)])), ...
                cat(3, ...
                    this.strandBreaks{i}([reshape(dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(1)), [], 1), 3*ones(nSites, 1)]),...
                    this.strandBreaks{i}([reshape(dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_RestrictionPositions(2)), [], 1), 4*ones(nSites, 1)]))];
            isPositionDamaged(isPositionStrandBreaks) = 0;
            isSiteUndamaged = ~any(any(isPositionDamaged, 2), 3);
            
            isSiteUnmethylated = ~any(isPositionMethylated, 2);
            isSitePolymerized = this.isRegionPolymerized(i, restrictionPosStrnds, 1, false, false, false);
            
            value = CircularSparseMat(...
                restrictionPosStrnds([isSiteUnmethylated; isSiteUnmethylated] & [isSiteUndamaged; isSiteUndamaged] & isSitePolymerized, :), ...
                true, [this.sequenceLen(i) this.nCompartments], 1);
        end
        
        function value = calcHemiunmethylatedMunIRMSites(this,i)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            dr = this.dnaRepair;
            nSites = size(dr.RM_MunI_RecognitionSites{i}, 1);
            methylationPosStrnds = [
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(1))   ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(1)) 3*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(2)) 2*ones(nSites, 1);
                dr.RM_MunI_RecognitionSites{i}(:, dr.RM_MunI_MethylatedPositions(2)) 4*ones(nSites, 1)];
            
            isPositionMethylated = reshape(this.damagedBases{i}(methylationPosStrnds), [], 2) == this.metabolite.m6ADIndexs;
            
            isSiteMethylated = all(isPositionMethylated, 2);
            isSiteUnmethylated = ~any(isPositionMethylated, 2);
            isSiteHemimethylated = ~isSiteMethylated & ~isSiteUnmethylated;
            
            value = CircularSparseMat(...
                methylationPosStrnds([isSiteHemimethylated; isSiteHemimethylated] & ~isPositionMethylated(:), :),  ...
                1, [this.sequenceLen(i) this.nCompartments], 1);
        end	
                
        function value = get.dryWeight(this)
			for i = 1:16
				if max(this.validated_polymerizedRegions(i), this.validated_damaged(i)) > this.validated_dryWeight(i)
					this.dryWeight{i} = this.calcDryWeight(i);
					this.validated_dryWeight(i) = this.validated(i);
				end
				value = this.dryWeight;
			end
				value = value{1} + value{2} + value{3} + value{4} + value{5} + value{6} + value{7} + value{8}...
					  + value{9} + value{10} + value{11} + value{12} + value{13} + value{14} + value{15} + value{16};
        end
        
        function value = calcDryWeight(this,i)
            import edu.jiangnan.fmme.util.ConstantUtil;

            %time
            numTime = size(this.abasicSites{i}, 3);
            
            %mass of undamaged DNA
            baseCounts = zeros(4, numTime);
            bonds = zeros(1, numTime);
			
            [positionsStrandTimes, lengths] = find(this.getStrandView(i, 'polymerizedRegions'));
			
            positionsStrandTimes = [positionsStrandTimes ones(size(positionsStrandTimes, 1), 3 - size(positionsStrandTimes, 2))];
            
            for j = 1:size(positionsStrandTimes, 1)
                if lengths(j) == this.sequenceLen(i)
                    baseCounts(:, positionsStrandTimes(j, 3)) = ...
                        baseCounts(:, positionsStrandTimes(j, 3)) + ...
                        getBaseCounts(this.sequence{i}, positionsStrandTimes(j, 2));
                else
                    baseCounts(:, positionsStrandTimes(j, 3)) = ...
                        baseCounts(:, positionsStrandTimes(j, 3)) + ...
                        this.sequence{i}.subsequenceBaseCounts(positionsStrandTimes(j,1) + (0:lengths(j)-1)', positionsStrandTimes(j, 2));
                end
                bonds(1, positionsStrandTimes(j,3)) = ...
                    bonds(1, positionsStrandTimes(j,3)) + lengths(j) - 1;
                
                if ...
                        lengths(j) == this.sequenceLen(i) || ...
                        (positionsStrandTimes(j, 1) + lengths(j) - 1 == this.sequenceLen(i) && ...
                        ismember([1 positionsStrandTimes(j, 2:end)], positionsStrandTimes, 'rows'))
                    bonds(1, positionsStrandTimes(j,3)) = ...
                        bonds(1, positionsStrandTimes(j,3)) + 1;
                end
            end
			
            value = this.metabolite.molecularWeights(this.metabolite.dnmpIndexs)' * baseCounts ...
                - (ConstantUtil.elements.H + ConstantUtil.elements.O) * bonds;
			            
            %mass represented by damage
            for k = 1:numTime
                %gap sites
                value(:,k) = value(:,k) - ...
                    (this.metabolite.molecularWeights(this.metabolite.dr5pIndexs) - ...
                    2 * ConstantUtil.elements.H) * collapse(this.gapSites{i}(:,:,k));
                
                %abasic sites
                [position, index] = find(this.abasicSites{i}(:,:,k));
                value(:,k) = value(:,k) - ...
                    length(index) * this.metabolite.molecularWeights(this.metabolite.waterIndexs) - ...
                    this.sequence{i}.subsequenceBaseCounts(position)' * ...
                    this.metabolite.molecularWeights(this.metabolite.dnmpIndexs);
                
                %damaged sugar-phosphates
                [~, index] = find(this.damagedSugarPhosphates{i}(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    length(index) * this.metabolite.molecularWeights(this.metabolite.dr5pIndexs);
                
                %damaged bases
                [position, index] = find(this.damagedBases{i}(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    this.sequence{i}.subsequenceBaseCounts(position)' * ...
                    this.metabolite.molecularWeights(this.metabolite.unmodifiedBaseIndexs);
                
                %intrastrand cross links
                [position, index] = find(this.intrastrandCrossLinks{i}(:,:,k));
                value(:,k) = value(:,k) + ...
                    sum(this.metabolite.molecularWeights(index)) - ...
                    this.sequence{i}.subsequenceBaseCounts([position(:,1) position(:,1) + 1], position(:,2))' * ...
                    this.metabolite.molecularWeights(this.metabolite.dnmpIndexs);
                
                %strand breaks
                value(:,k) = value(:,k) + ...
                    (ConstantUtil.elements.H + ConstantUtil.elements.O)*...
                    collapse(this.strandBreaks{i}(:,:,k));
            end
            
            value = permute(value, [1 3 2]);
            value = value / ConstantUtil.nAvogadro;			
        end
    end
    
    %helper methods
    methods (Static)
        function value = shiftCircularSparseMatBase3Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBase3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBase5Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBase5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBond3Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBond3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = unshiftCircularSparseMatBond3Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.unshiftPositionsStrandsBond3Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = shiftCircularSparseMatBond5Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.shiftPositionsStrandsBond5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function value = unshiftCircularSparseMatBond5Prime(sparseMat, varargin)
            import edu.jiangnan.fmme.cell.sim.state.Chromosome;
            import edu.jiangnan.fmme.util.CircularSparseMat;
            
            if ~isa(sparseMat,'CircularSparseMat') || ~isDimCircular(sparseMat, 1)
                throw(MException('Chromosome:error','sparseMat must be a circular sparse mat with a circular first dimension'));
            end
            
            [subs, vals] = find(sparseMat);
            subs = Chromosome.unshiftPositionsStrandsBond5Prime(subs, varargin{:});
            value = CircularSparseMat(subs, vals, size(sparseMat), 1);
        end
        
        function positionsStrands = shiftPositionsStrandsBase3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) + lengths;
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) - lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBase5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) - lengths;
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) + lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBond3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) - lengths;
        end
        
        function positionsStrands = unshiftPositionsStrandsBond3Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1; 
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==0, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==0, 1) + lengths;
        end
        
        function positionsStrands = shiftPositionsStrandsBond5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1;
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) - lengths;
        end
        
        function positionsStrands = unshiftPositionsStrandsBond5Prime(positionsStrands, lengths)
            if nargin < 2
                lengths = 1;
            end
            
            positionsStrands(mod(positionsStrands(:,2),2)==1, 1) = ...
                positionsStrands(mod(positionsStrands(:,2),2)==1, 1) + lengths;
        end
        
        function [footprint3Prime, footprint5Prime] = calculateFootprintOverhangs(footprint)
            footprint5Prime = ceil((footprint-1)/2);
            footprint3Prime = footprint - 1 - footprint5Prime;
        end
    end
    
    %printing
    methods
        %print state
        function disp(this)
            %superclass method
            this.disp@edu.jiangnan.fmme.cell.sim.CellState();
            %numbers of DNA damages
			for i = 1:16
            fprintf('%24s\t%4s\n','Damage Type','No.');
            fprintf('%24s\t%4s\n',repmat('=',1,24), repmat('=',1,4));
            fprintf('%24s\t%4d\n','Gap Sites',               collapse(this.gapSites{i}));
            fprintf('%24s\t%4d\n','Abasic Sites',            collapse(this.abasicSites{i}));
            fprintf('%24s\t%4d\n','Damaged sugar phosphates',collapse(this.damagedSugarPhosphates{i}));
            fprintf('%24s\t%4d\n','Damaged bases',           collapse(this.damagedBases{i}));
            fprintf('%24s\t%4d\n','Intrastrand cross links', collapse(this.intrastrandCrossLinks{i}));
            fprintf('%24s\t%4d\n','Strand breaks',           collapse(this.strandBreaks{i}));
            fprintf('%24s\t%4d\n','Holliday junctions',      collapse(this.hollidayJunctions{i}));
            fprintf('\n');
           
            %strand break classification used in track structure models
			sbc = this.strandBreakClassification{1}();
            fprintf('%24s\t%4s\n','Strand Break','No.')
            fprintf('%24s\t%4s\n',repmat('=',1,24), repmat('=',1,4))
            fprintf('%24s\t%4d\n','NB',    sbc(this.strandBreakClassification_index_NB));
            fprintf('%24s\t%4d\n','SSB',   sbc(this.strandBreakClassification_index_SSB));
            fprintf('%24s\t%4d\n','SSB+',  sbc(this.strandBreakClassification_index_SSB_));
            fprintf('%24s\t%4d\n','2SSB',  sbc(this.strandBreakClassification_index_2SSB));
            fprintf('%24s\t%4d\n','DSB',   sbc(this.strandBreakClassification_index_DSB));
            fprintf('%24s\t%4d\n','DSB+',  sbc(this.strandBreakClassification_index_DSB_));
            fprintf('%24s\t%4d\n','DSB++', sbc(this.strandBreakClassification_index_DSB__));
            fprintf('\n');
            
            %list of damaged genes
            fprintf('%13s\n','Damaged Genes');
            fprintf('%13s\n',repmat('=',1,13));
            damagedGeneWholeCellModelIDs{i} = this.gene.wholeCellModelIDs(~all(this.accessibleGenes{i},2));
            damagedGeneWholeCellNames{i}    = this.gene.names(~all(this.accessibleGenes{i},2));
            for j=1:length(damagedGeneWholeCellModelIDs{i});
                damagedGeneWholeCellModelID{i} = damagedGeneWholeCellModelIDs{i}{j};
                damagedGeneWholeCellName{i}    = damagedGeneWholeCellNames{i}{j};
                if length(damagedGeneWholeCellModelID{i})>12; damagedGeneWholeCellModelID{i}=[damagedGeneWholeCellModelID{i}(1:8) ' ...']; end;
                if length(damagedGeneWholeCellName{i})>32; damagedGeneWholeCellName{i}=[damagedGeneWholeCellName{i}(1:28) ' ...']; end;
                fprintf('%12s\t%32s\n',damagedGeneWholeCellModelID{i},damagedGeneWholeCellName{i});
            end
            fprintf('\n');
            
            %list of damaged transcription units
            fprintf('%27s\n','Damaged Transcription Units');
            fprintf('%27s\n',repmat('=',1,27));
            damagedTranscriptionUnitWholeCellModelIDs{i} = this.transcriptionUnitWholeCellModelIDs{i}(~all(this.accessibleTranscriptionUnits{i},2));
            damagedTranscriptionUnitWholeCellNames{i}    = this.transcriptionUnitNames{i}(~all(this.accessibleTranscriptionUnits{i},2));
            for j=1:length(damagedTranscriptionUnitWholeCellModelIDs{i});
                damagedTranscriptionUnitWholeCellModelID{i} = damagedTranscriptionUnitWholeCellModelIDs{i}{j};
                damagedTranscriptionUnitWholeCellName{i}    = damagedTranscriptionUnitWholeCellNames{i}{j};
                if length(damagedTranscriptionUnitWholeCellModelID{i})>12; damagedTranscriptionUnitWholeCellModelID{i}=[damagedTranscriptionUnitWholeCellModelID{i}(1:8) ' ...']; end;
                if length(damagedTranscriptionUnitWholeCellName{i})>32; damagedTranscriptionUnitWholeCellName{i}=[damagedTranscriptionUnitWholeCellName{i}(1:28) ' ...']; end;
                fprintf('%12s\t%32s\n',damagedTranscriptionUnitWholeCellModelID{i},damagedTranscriptionUnitWholeCellName{i});
            end
            fprintf('\n');
			end
        end
    end
end
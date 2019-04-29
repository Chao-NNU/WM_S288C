%DNASupercoiling
%
% @wholeCellModelID Process_DNASupercoiling
% @name             DNA Supercoiling
% @description
%   Biology
%   ===========
%   DNA gyrase and topoisomerase IV each use 2 ATP to induce 2
%   negative
%   supercoils each time they act. Occasionally, topoisomerase IV can also
%   induce positive supercoils but that is not considered in our model.
%   Topoisomerase I acts to induce positive supercoils. Gyrase, topoisomerase
%   IV, and topoisomerase I act at rates of 1.2, 2.5, and 1 strand passing
%   events per second respectively.
% 
%   We use a calculation of the DNA's linking number (LK) in order to track
%   the supercoiling of the DNA.  The &delta;LK is the difference between the
%   current level of DNA supercoiling and the relaxed level of DNA
%   supercoiling. The LKrelaxed is defined as number of base pairs/10.5, in which
%   10.5 is the number of bases per turn in a relaxed double helix. As the
%   replication loops move, the LKcurrent deviates from this relaxed state,
%   and gyrases and topoisomerases help bring the DNA back to the relaxed
%   state. The superhelical density, or specific linking number density,
%   &sigma;<sub>sp</sub>, is defined as the:
%      (LKcurrent - LKrelaxed)/LKrelaxed. 
%   The activity of gyrases and topoisomerases depends on the &sigma;<sub>sp</sub>
%   of the DNA. Although there is likely a more complex relationship between the
%   DNA supercoiling and enzyme activity, we model the enzyme activity as a
%   combination of step functions and logistic functions. 
%   Topoisomerase IV can only act if &sigma;<sub>sp</sub> is
%   higher than 0. Topoisomerase I can only act if &sigma;<sub>sp</sub> is lower
%   than zero. Gyrase can only act if &sigma;<sub>sp</sub> is higher than -0.1.
%   Within, the regions where gyrase and topoI are allowed to act, we apply
%   a logistic function describing the probability of activity, centered
%   around equilibrium sigma, -0.06. For TopoI, the probability of activity
%   approaches 1 as the sigma gets more positive. For gyrase, the
%   probability of activity approaches 1 as
%   
%   We model up to three regions on the chromosomes and track their LKs
%   separately. Before replication, "unreplicated DNA" is the only region
%   present. The enzymes may act on the chromosome affecting the LKcurrent,
%   and experimentally it has been shown that the enzymes would obtain an
%   LKcurrent such that the steady state &sigma;<sub>sp</sub> is =-0.06. As the
%   replication loop progresses, "unreplicated DNA" is the region downstream of
%   the two replication loops.  The number of bases in this region decreases
%   during replication, meaning that the LKrelaxed decreases. The LKcurrent then,
%   is too high, and must be brought down towards the LKrelaxed by inducing
%   negative supercoils.
%   
%   The second and third regions are the replicated DNA upstream of the
%   replication loops on each of the two chromosomes. As new uncoiled DNA is
%   formed, it is already coiled, but gyrases and topoisomerases can continue
%   to act on this DNA, ideally maintaining the steady state &sigma;<sub>sp</sub>.
%   After replication is complete, these two regions are the only two that exist.
%   
%   It is essential that the &detlta;LK in the region downstream of the replication
%   loops be brought down to 0 by the end of replication.
%   
%   Another consideration is the processivity of the enzymes when bound to the
%   DNA. Topoisomerase IV is highly processive and will stay bound to the DNA
%   as long as the &sigma;<sub>sp</sub> is greater than zero. Topoisomerase I is
%   not highly processive, and essentially acts on the DNA and falls back off
%   right away gyrase will stay bound to the DNA for about 30-60 seconds. We model
%   gyrase processivity as a poisson distribution with &lambda; = 45 seconds.
%   
%   In this process, we go through all free gyrases and topoisomerases in
%   random order,  determine what regions they can bind in, and randomly bind
%   them to a large enough open position on the DNA. We adjust the linking
%   number based on all enzyme actions, and account for the usage of ATP. We
%   also track the processivity of gyrase and topoisomerase IV, and unbind
%   them from the DNA when appropriate.
%   
%   If replication is in progress, we first knock off any enzymes that the
%   replication loop collides into.
%
%   There is an effect of supercoiling on the probabilities of gene 
%   transcription (Peter 2004). While fold change at differnt sigmas have
%   been calculated for many E. coli genes, here, for simplicity, we are
%   only including the effects on the 5 supercoiling genes: gyrB, gyrA,
%   parC, parE, and topA. These genes exist in 3 transcription units. Peter
%   2004 has data for the fold change of expression of each of these genes
%   are various values of sigma (ranging from sigma -0.06-0.02) at various 
%   experimental conditions. Due to the limited data and large variation within
%   the data, we have decided to use a linear fit of the fold change data and 
%   extrapolate the linear fit within the sigma range of -0.08 to 0.07 where 
%   the linear fit seems reasonable. Outside of this range, we estimate a
%   constant fold change of expression. While there is a separate set of
%   data for each of the 5 genes, our model requires a single probability
%   of transcription for each transcription unit. The data for the first
%   gene in each transcription unit is used (gyrB, and parE). TopA is
%   transcribed with genes that are not supercoiling related. The
%   expression of those genes (MG_119, 120, 121) will also get affected by 
%   this process. The general trend is that gyrase and topoIV will have an
%   increased expression when sigma is higher than the equilibrium, and
%   that topoI will have a higher expression when sigma is lower than the
%   equilibrium. Fold changes for gyrase will vary between 0.14 and 6.57.
%   Fold changes for topoIV will vary between 0.98 and 1.14. Fold changes
%   for topoI will vary between 0.042 and 1.15. 
%   
%   References
%   ===============
%   1. Ullsperger, C., Cozzarelli, N.R. (1996). Contrasting enzymatic activities
%      of topoisomerase IV and DnA gyrase from Escherichia coli. Journal of Bio
%      Chem 271: 31549-31555. [PUB_0236]
%   2. Dekker, N.H., Viard, T., Bouthier de la Tour, C., Duguet, M., Bensimon,
%      D., Croquette, V. (2003). Thermophilic Topoisomerase I on a single DNA
%      molecule. Journal of molecular biology 329: 271-282. [PUB_0502]
%   3. Gore, J., Bryant, Z., Stone, M.D., Nollmann, M., Cozzarelli, N.R.,
%      Bustamante, C. (2006). Mechanochemical analysis of DNA gyrase using rotor
%      bead tracking. Nature 439: 100-104. [PUB_0751]
%   4. Bates, A. (2006). DNA Topoisomerases: Single Gyrase Caught in the Act.
%      Current Biology 16: 204-206. [PUB_0752]
%   5. Peng, H., Marians, K.J. (1995). The Interaction of Escherichia coli
%      Topoisomerase IV with DNA. Journal of Biological Chemistry 42:
%      25286?25290. [PUB_0694]
%   6. Wang, J. (1996) DNA Topoisomerases. Annual Reviews 65: 635-692. [PUB_0693]
%   7. Peter, B.J., Arsuaga, J., Breier, A.M., Khodursky, A.B., Brown,
%      P.O., Cozzarelli, N.R. (2004) Genomic transcriptional response to loss
%      of chromosomal supercoiling in Escherichia coli. Genome Biology 5:
%      1-13. [PUB_0920]
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 11/18/2010
classdef DNASupercoiling < edu.jiangnan.fmme.cell.sim.Process & edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect
    %property annotations
    properties (Constant)
        optionNames__              = {}; %names of option properties
        fixedConstantNames__       = {   %names of fixed constant properties
            'gyraseMeanDwellTime';
            'gyraseSigmaLimit';
            'topoISigmaLimit';
            'topoIVSigmaLimit';
            'gyraseDeltaLK';
            'topoIDeltaLK';
            'topoIVDeltaLK';
            'gyraseActivityRate';
            'topoIActivityRate';
            'topoIVActivityRate';
            'gyraseATPCost';
            'topoIATPCost';
            'topoIVATPCost';
            'foldChangeSlopes'; 
            'foldChangeIntercepts';
            'foldChangeLowerSigmaLimit';
            'foldChangeUpperSigmaLimit'; 
            'numTranscriptionUnits';
            'tuIndexs';
            'tuCoordinates';
            };
        fittedConstantNames__      = {}; %names of fitted constant properties
        localStateNames__          = {}; %names of simulation state properties redundant with timecourses in this or other processes or the simulation
    end

    %IDs, names, and local indices
    properties
        stimuliWholeCellModelIDs = {};   %whole cell model IDs of stimuli
        substrateWholeCellModelIDs = {   %whole cell model IDs of substrates
            'ATP'; 'ADP'; 'PI'; 'H2O'; 'H'};
        substrateIndexs_atp       = 1;
        substrateIndexs_adp       = 2;
        substrateIndexs_phosphate = 3;
        substrateIndexs_water     = 4;
        substrateIndexs_hydrogen  = 5;

        enzymeWholeCellModelIDs = {      			%whole cell model IDs of enzymes
            'YPL024W_YMR190C_YLR234W_TRIMER';       %DNA gyrase
            'YNL088W_MONOMER';      				%DNA topoisomerase IV
            'YOL006C_MONOMER'};           			%DNA topoisomerase I
        enzymeIndexs_gyrase = 1;
        enzymeIndexs_topoIV = 2;
        enzymeIndexs_topoI  = 3;
    end
    
    %fixed biological constants
    properties
        gyraseMeanDwellTime       %Mean dwell time of gyrase to bind DNA [45s; PUB_0753]
        
        gyraseSigmaLimit          %superhelical density (\sigma_{sp}) above which gyrase acts [-0.1; PUB_0753]
        topoISigmaLimit           %superhelical density (\sigma_{sp}) above which topoisomerase I acts [0; PUB_0502]
        topoIVSigmaLimit          %superhelical density (\sigma_{sp}) below which topoisomerase IV acts [0; PUB_0753]
        
        gyraseDeltaLK             %number of supercoils introduced by each gyrase strand passing event [-2; PUB_0695]
        topoIDeltaLK              %number of supercoils introduced by each topoisomerase I strand passing event [1; PUB_0502]
        topoIVDeltaLK             %number of supercoils introduced by each topoisomerase IV strand passing event [-2]
        
        gyraseActivityRate        %strand passing events per second per enzyme [1.2; PUB_0753]
        topoIActivityRate         %strand passing events per second per enzyme [1.0; PUB_0502]
        topoIVActivityRate        %strand passing events per second per enzyme [2.5; PUB_0750]
        
        gyraseATPCost             %ATP hydrolyzed per gyrase-catalyzed strand passing event [2]
        topoIATPCost              %ATP hydrolyzed per topoisomerase I-catalyzed strand passing event [0]
        topoIVATPCost             %ATP hydrolyzed per topoisomerase IV-catalyzed strand passing event [2]
        
        topoILogisiticConst       %fittable parameter controlling the steapness of the the logistic function of topoI activity around equilibrium sigma [100]
        gyrLogisiticConst         %fittable parameter controlling the steapness of the the logistic function of gyrase activity around equilibrium sigma [-100]
        
        enzymeProperties          %struct built by buildEnzymeProperties
        
        foldChangeSlopes  		  %slopes of linear fit of sigma-gene expression data [gyrase, topo IV, topo I] [PUB_0920] 
        foldChangeIntercepts      %y-intercepts of linear fit of sigma-gene expression data [gyrase, topo IV, topo I] [PUB_0920]
        foldChangeLowerSigmaLimit %lower bound of sigma for applying linear fit (fittable parameter)
        foldChangeUpperSigmaLimit %upper bound of sigma for applying linear fit (fittable parameter)
        
        numTranscriptionUnits     %number of transcription units
        tuIndexs 			      %transcription units whose transcription probs will be affected by supercoiling [gyrase, topo IV, topo I]
        tuCoordinates             %start coords of transcription units whose transcription probs will be affected by supercoiling [gyrase, topo IV, topo I]
    end
    
    %references to cell state
    properties
        rnaPolymerase             %reference to RNA polymerase state
    end

    %constructor
    methods
        function this = DNASupercoiling(wholeCellModelID, name)
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
        function initializeConstants(this, knowledgeBase, simulation, varargin)
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.Process(...
                knowledgeBase, simulation, varargin{:});
            this.initializeConstants@edu.jiangnan.fmme.cell.sim.ChromosomeProcessAspect(...
                knowledgeBase, simulation, varargin{:});            
            
            %c = this.chromosome;
			c = knowledgeBase.transcriptionUnits;
			
            this.enzymeProperties = this.buildEnzymeProperties();
            
            %get transcription unit indices			
            
			[~, geneIndexs] = ismember({'YOL006C', 'YPL024W', 'YNL088W'}, this.gene.wholeCellModelIDs);
            [this.tuIndexs, ~] = find(this.rna.nascentRNAGeneComposition(geneIndexs, :)');   
            %this.numTranscriptionUnits = numel(c.transcriptionUnitStartCoordinates);
			this.numTranscriptionUnits = numel(c);
			s = [c.startCoordinate]';
            %this.tuCoordinates = c.transcriptionUnitStartCoordinates(this.tuIndexs);
			this.tuCoordinates = s(this.tuIndexs);
        end
    end

    %model
    methods        
        %Calculate sufficient metabolites and enzymes to achieve the observed
        %mean superhelical density after replication of the daughter chromosome:
        %
        %  (LK after replication) = 2 * (LK before replication)
        %                         = (LK before replication)
        %                            + (links introduced by DNA pol)
        %                            - (links pushed ahead of replisome and resolved by gyrase)
        %                            - (links introduced behind replisome by gryase to return to steady-state superhelical density)
        %               2*LK_{ss} = LK_{ss} + 2*LK_{0} - LK_{ss} - 2*(LK_{0} - LK_{ss})
        %
        %That is calculate:
        %- contribution to FBA objective
        %- minimum expression consistent with cell cycle length
        %
        %FitConstants accounts for topoisomerase I / gyrase expression balance
        %required to achieve the observed steady state superhelical density.
        function [bmProd, byProd, minEnzExp, maxEnzExp] = calcResourceRequirements_LifeCycle(this, constants, states)
            %initialize
            bmProd = zeros(size(this.substrateWholeCellModelIDs));
            byProd = zeros(size(this.substrateWholeCellModelIDs));
            minEnzExp = zeros(size(this.enzymeWholeCellModelIDs));
            maxEnzExp = Inf(size(this.enzymeWholeCellModelIDs));
            
            %references
            c = this.chromosome;
            t = constants.states.Time;
            
            %number of links to be resolved by gyrase, and number of gyrase
            %link resolution reactions
            %- links that need to be resolved because of replication
            %- links that need to be resolved to counter-balance the effect
            %  of topoisomerase I
            
            for i = 1:16
				nRxn(i) = 0;
				LK_0(i)  =  [c.sequenceLen(i) / c.relaxedBasesPerTurn(i)];
				%LK_0  = LK_0(:,4);
				LK_ss(i) = (1 + c.equilibriumSuperhelicalDensity) * LK_0(i);
				dLK(i) = - LK_ss(i) - 2*(LK_0(i) - LK_ss(i));
				nRxn(i) = nRxn(i) + dLK(i) / this.gyraseDeltaLK;
			end
			nRxn = nRxn(1) + nRxn(2) + nRxn(3) + nRxn(4) + nRxn(5) + nRxn(6) + nRxn(7) + nRxn(8)...
				 + nRxn(9) + nRxn(10) + nRxn(11) + nRxn(12) + nRxn(13) + nRxn(14) + nRxn(15) + nRxn(16);
            
            topoIProps = this.enzymeProperties(this.enzymeIndexs_topoI);
            gyrProps = this.enzymeProperties(this.enzymeIndexs_gyrase);
            nRxn = nRxn + 0.25 * ... %factor of 0.25 added based on observing simulations
                abs(states.monomers(this.enzymeGlobalIndexs(this.enzymeIndexs_topoI)) * ...
                topoIProps.activityRate * topoIProps.deltaLK / gyrProps.deltaLK);
            
            %substrates and byproducts
            nATP = nRxn * this.enzymeProperties(this.enzymeIndexs_gyrase).atpCost;
			nATP = sum(nATP);
            bmProd(this.substrateIndexs_atp)       = nATP;
            bmProd(this.substrateIndexs_water)     = nATP;
            byProd(this.substrateIndexs_adp)       = nATP;
            byProd(this.substrateIndexs_phosphate) = nATP;
            byProd(this.substrateIndexs_hydrogen)  = nATP;
            
            %Gyrase
            minEnzExp(this.enzymeIndexs_gyrase) = ...
                sum(nRxn / t.replicationDuration / this.gyraseActivityRate / ...
                exp(log(2) * t.replicationInitiationDuration / t.cellCycleLength));
        end

        %initialization:
        %- Simulation initializeState method initializes 1
        %  chromosome with linkingNumber such that the superhelical density is
        %  equal to the equilibriumSuperhelicalDensity
        %- Here we bind all gyrase to the  
        function initializeState(this)
            %bind gyrase
			c = this.chromosome;
			for i= 1:16
            nGyraseBinding = this.randStream.stochasticRound(...
                this.enzymes(this.enzymeIndexs_gyrase) * ...
                (1 - 1/this.gyraseMeanDwellTime/this.stepSizeSec));
				
            this.bindProteinToChromosomeStochastically(i, this.enzymeIndexs_gyrase, ...
                nGyraseBinding);
            
            %Calculate effects of supercoiling on probability of transcription
			this.rnaPolymerase.supercoilingBindingProbFoldChange{i} = this.calcRNAPolymeraseBindingProbFoldChange(i);
			end
        end
        
        %resource requirements
        function result = calcResourceRequirements_Current(this)
            %chromosome reference
            c = this.chromosome;
			for i = 1:16
				%initialize
				result = zeros(size(this.substrates));
            
				%ATP needed for gyrase activity
				[~, sigmas{i}] = find(c.superhelicalDensity{i});
				if isempty(sigmas{i})
					sigmas{i} = 0;
				end
            end
            enzProps = this.enzymeProperties(this.enzymeIndexs_gyrase);
            nATP = enzProps.atpCost * ceil(...
                (this.enzymes(this.enzymeIndexs_gyrase) + this.boundEnzymes(this.enzymeIndexs_gyrase)) * ...
                enzProps.activityRate *  max(enzProps.probOfActivityFunc(...
				(sigmas{1} + sigmas{2} + sigmas{3} + sigmas{4} + sigmas{5} + sigmas{6} + sigmas{7} + sigmas{8}...
				 + sigmas{9} + sigmas{10} + sigmas{11} + sigmas{12} + sigmas{13} + sigmas{14} + sigmas{15} + sigmas{16}), c)));
            
            result(this.substrateIndexs_atp)   = nATP;
            result(this.substrateIndexs_water) = nATP;
        end

        %simulation
        function evolveState(this)
            import edu.jiangnan.fmme.util.CircularSparseMat;
            import edu.jiangnan.fmme.util.ComputationUtil;
            
            c = this.chromosome;
			
            for i = 1: 16
			
            %get enzyme properties
            enzProps = this.enzymeProperties;

            % Calculate sigmas
            [tmpPosStrands{i}, lengths{i}] = find(c.doubleStrandedRegions{i});
            tmpIdxs{i} = find(mod(tmpPosStrands{i}(:, 2), 2));
            dsPosStrands{i} = tmpPosStrands{i}(tmpIdxs{i}, :);
            lengths{i} = lengths{i}(tmpIdxs{i}, 1);
            
            [~, linkingNumbers{i}] = find(c.linkingNumbers{i});
            linkingNumbers{i} = linkingNumbers{i}(tmpIdxs{i}, 1);
			
            sigmas{i} = (linkingNumbers{i} - lengths{i}/c.relaxedBasesPerTurn(i)) ./ (lengths{i}/c.relaxedBasesPerTurn(i));
			
            % Determine which regions accommodate which proteins.
            legal{i} = [...
                sigmas{i} > this.gyraseSigmaLimit ...
                sigmas{i} > this.topoIVSigmaLimit ...
                sigmas{i} < this.topoISigmaLimit];

				% Unbind topoIV where it can't act.
            this.releaseProteinFromChromosome(i,...
                this.enzymeIndexs_topoIV, Inf,...
                dsPosStrands{i}(legal{i}(:, this.enzymeIndexs_topoIV), :),... %protected regions
                lengths{i}(legal{i}(:, this.enzymeIndexs_topoIV)));
            
            % Probabilistically unbind gyrase.
            this.releaseProteinFromChromosome(i, this.enzymeIndexs_gyrase, 1 / this.gyraseMeanDwellTime, [], []);
            
            % Bind enzymes in regions with accommodating sigmas.
            order = this.randStream.randperm(length(this.enzymes));
			
            transientlyBinding{i} = zeros(size(legal{i}));
			
			%keyboard
			
			for m = 1:length(order)
                j = order(m);
                
                if enzProps(j).meanDwellTime == 0
                    if numel(lengths{i}) == 1 && lengths{i} == c.sequenceLen(i)
                        [~, mons{i}] = find(c.monomerBoundSites{i});
                        [~, cpxs{i}] = find(c.complexBoundSites{i});
                        space{i} = lengths{i} ...
                            - sum(c.monomerDNAFootprints(mons{i}, 1)) ...
                            - sum(c.complexDNAFootprints(cpxs{i}, 1)) ...
                            - nnz(c.damagedSites{i});
                        transientlyBinding{i}(1, j) = min(this.enzymes(j, 1), space{i} / this.enzymeDNAFootprints(j, 1));
                    elseif numel(lengths{i}) == 2 && all(lengths{i} == c.sequenceLen(i))
                        [monPosStrnds{i}, mons{i}] = find(c.monomerBoundSites{i});
                        [cpxPosStrnds{i}, cpxs{i}] = find(c.complexBoundSites{i});
                        space{i} = zeros(2, 1);
                        space{i}(1) = lengths{i}(1) ...
                            - sum(c.monomerDNAFootprints(mons{i}(monPosStrnds{i}(:, 2) <=2 ), 1)) ...
                            - sum(c.complexDNAFootprints(cpxs{i}(cpxPosStrnds{i}(:, 2) <=2 ), 1)) ...
                            - nnz(c.damagedSites{i}(:, 1:2));
                        space{i}(2) = lengths{i}(2) ...
                            - sum(c.monomerDNAFootprints(mons{i}(monPosStrnds{i}(:, 2) >2 ), 1)) ...
                            - sum(c.complexDNAFootprints(cpxs{i}(cpxPosStrnds{i}(:, 2) >2 ), 1)) ...
                            - nnz(c.damagedSites{i}(:, 3:4));
                        transientlyBinding{i}(:, j) = min(this.enzymes(j, 1), sum(space{i}) / this.enzymeDNAFootprints(j, 1)) * ...
                            space{i} / sum(space{i});
                        if this.randStream.rand < 0.5
                            transientlyBinding{i}(1, j) = ComputationUtil.roundHalfUp(transientlyBinding{i}(1, j));
                            transientlyBinding{i}(2, j) = ComputationUtil.roundHalfDown(transientlyBinding{i}(2, j));
                        else
                            transientlyBinding{i}(2, j) = ComputationUtil.roundHalfUp(transientlyBinding{i}(2, j));
                            transientlyBinding{i}(1, j) = ComputationUtil.roundHalfDown(transientlyBinding{i}(1, j));
                        end
                    else
                        [monPosStrnds{i}, mons{i}] = find(c.monomerBoundSites{i});
                        [cmpPosStrnds{i}, cmps{i}] = find(c.complexBoundSites{i});
                        dmgPosStrnds{i} = find(c.damagedSites{i});
                        monPos{i} = monPosStrnds{i}(:, 1);
                        cmpPos{i} = cmpPosStrnds{i}(:, 1);
                        dmgPos{i} = dmgPosStrnds{i}(:, 1);
                        monChrs{i} = ceil(monPosStrnds{i}(:, 2) / 2);
                        cmpChrs{i} = ceil(cmpPosStrnds{i}(:, 2) / 2);
                        dmgChrs{i} = ceil(dmgPosStrnds{i}(:, 2) / 2);
                        
                        space{i} = zeros(size(lengths{i}));
                        for k = 1:numel(lengths{i})
                            if ~legal{i}(k, j)
                                continue;
                            end
                            
                            monTfs{i} = monPos{i} >= dsPosStrands{i}(k, 1) & monPos{i} <= dsPosStrands{i}(k, 1) + lengths{i}(k, 1) - 1;
                            cmpTfs{i} = cmpPos{i} >= dsPosStrands{i}(k, 1) & cmpPos{i} <= dsPosStrands{i}(k, 1) + lengths{i}(k, 1) - 1;
                            dmgTfs{i} = dmgPos{i} >= dsPosStrands{i}(k, 1) & dmgPos{i} <= dsPosStrands{i}(k, 1) + lengths{i}(k, 1) - 1;
                            monTfs{i}(monTfs{i}) = monChrs{i}(monTfs{i}, 1) == dsPosStrands{i}(k, 2);
                            cmpTfs{i}(cmpTfs{i}) = cmpChrs{i}(cmpTfs{i}, 1) == dsPosStrands{i}(k, 2);
                            dmgTfs{i}(dmgTfs{i}) = dmgChrs{i}(dmgTfs{i}, 1) == dsPosStrands{i}(k, 2);
                            
                            space{i}(k, 1) = ...
                                + lengths{i}(k, 1) ...
                                - sum(c.monomerDNAFootprints(mons{i}(monTfs{i}, 1), 1)) ...
                                - sum(c.monomerDNAFootprints(cmps{i}(cmpTfs{i}, 1), 1)) ...
                                - sum(dmgTfs{i});
                        end
                        
                        if numel(lengths{i}) == 1
                            transientlyBinding{i}(:, j) = min(this.enzymes(j, 1), space{i} / this.enzymeDNAFootprints(j, 1));
                        else
                            transientlyBinding{i}(:, j) = space{i} / sum(space{i}) * min(this.enzymes(j, 1), sum(space{i}) / this.enzymeDNAFootprints(j, 1));
                        end
                    end
                elseif any(legal{i}(:, j))
                    this.bindProteinToChromosomeStochastically(...
                       i, j, [], dsPosStrands{i}(legal{i}(:,j),:), lengths{i}(legal{i}(:,j)));
                end
            end

            % Calculate new linking numbers and do substrate accounting.
            enzProps = enzProps(this.randStream.randperm(length(enzProps)));
            for m = 1:size(dsPosStrands{i}, 1)
                len = lengths{i}(m);
                for j = 1:length(enzProps)
                    if legal{i}(m, enzProps(j).idx)
                        if enzProps(j).meanDwellTime == 0
                            nBound{i} = transientlyBinding{i}(m, enzProps(j).idx);
                        else
                            nBound{i} = size(this.findProteinInRegion(...
                                i, dsPosStrands{i}(m, 1), dsPosStrands{i}(m, 2), len, enzProps(j).idx), 1);
                        end
                        % use logistic function to determine the
                        % probability of enzyme activity
                        % Use sigmas from beginning of timestep, before any
                        % enzymes have acted
                        probOfActivity{i} = enzProps(j).probOfActivityFunc(sigmas{i}(m),c);
                        nStrandPassingEvents{i} = min([
                            this.randStream.stochasticRound(nBound{i} * enzProps(j).activityRate *  probOfActivity{i});
                            fix(this.substrates([this.substrateIndexs_atp; this.substrateIndexs_water]) / enzProps(j).atpCost)
                            ]);
                        %adjust linking number
                        linkingNumbers{i}(m) = linkingNumbers{i}(m) + enzProps(j).deltaLK * nStrandPassingEvents{i};
                        
                        %substrate accounting
                        nATP{i} = nStrandPassingEvents{i} * enzProps(j).atpCost;
						%{
						nATP = nATP{1} + nATP{2} + nATP{3} + nATP{4} + nATP{5} + nATP{6} + nATP{7} + nATP{8} + nATP{9} + nATP{10}...
							    + nATP{11} + nATP{12} + nATP{13} + nATP{14} + nATP{15} + nATP{16};
						%}		
                        this.substrates(this.substrateIndexs_atp)       = this.substrates(this.substrateIndexs_atp)       - nATP{i};
                        this.substrates(this.substrateIndexs_water)     = this.substrates(this.substrateIndexs_water)     - nATP{i};
                        this.substrates(this.substrateIndexs_adp)       = this.substrates(this.substrateIndexs_adp)       + nATP{i};
                        this.substrates(this.substrateIndexs_phosphate) = this.substrates(this.substrateIndexs_phosphate) + nATP{i};
                        this.substrates(this.substrateIndexs_hydrogen)  = this.substrates(this.substrateIndexs_hydrogen)  + nATP{i};
                    end
                end
            end
            c.linkingNumbers{i} = CircularSparseMat(...
                [dsPosStrands{i}; dsPosStrands{i}(:, 1) dsPosStrands{i}(:, 2)+1], [linkingNumbers{i}; linkingNumbers{i}], ...
                [c.sequenceLen(i) c.nCompartments], 1); 
			
            this.rnaPolymerase.supercoilingBindingProbFoldChange{i} = this.calcRNAPolymeraseBindingProbFoldChange(i, dsPosStrands{i},...
																	  lengths{i}, linkingNumbers{i}, c.transcriptionUnitStartCoordinates{i},...
																	  c.transcriptionUnitStartCoordinates{i});
			end
		end
        
        function foldChange = calcRNAPolymeraseBindingProbFoldChange(this, i, dsPosStrands, lengths, linkingNumbers,...
							  transcriptionUnitStartCoordinates, relaxedBasesPerTurn)
            c = this.chromosome;
            if nargin == 2
                [tmpPosStrands, lengths] = find(c.doubleStrandedRegions{i});
                tmpIdxs = find(mod(tmpPosStrands(:, 2),2));
                dsPosStrands = tmpPosStrands(tmpIdxs, :);
                lengths = lengths(tmpIdxs, 1);
                [~, linkingNumbers] = find(c.linkingNumbers{i});
                linkingNumbers = linkingNumbers(tmpIdxs, 1);
			end
			
            foldChange = ones(size(c.transcriptionUnitStartCoordinates{i}, 1), 2); %initialize fold change
            sigmas = (linkingNumbers - lengths/c.relaxedBasesPerTurn(i)) ./ (lengths/c.relaxedBasesPerTurn(i));
			
			this.tuIndexs= [6, 24, 87];
			
			if i == 16
				for k = 2
				%for k = 1:length(this.tuCoordinates)% change from i to k
					%find double stranded region where transcription unit lies
					rgnIdxs = find(...
						dsPosStrands(:,1)               <= this.tuCoordinates(k) & ...
						dsPosStrands(:,1) + lengths - 1 >= this.tuCoordinates(k));
					%rgnIdxs = 1;
					for j = 1:numel(rgnIdxs)
						%piecewise linear function of sigma of double stranded region
						%where transcription unit lies
					
						thresholdedSigma = min(...
							max(sigmas(rgnIdxs(j)), this.foldChangeLowerSigmaLimit),...
							this.foldChangeUpperSigmaLimit);
					
						foldChange(...
							this.tuIndexs(k), ceil(dsPosStrands(rgnIdxs(j), 2)/2)) = ...
							this.foldChangeIntercepts(k) + this.foldChangeSlopes(k) * thresholdedSigma;
					end
				end
			elseif i == 15
				for k = 1
				rgnIdxs = find(...
						dsPosStrands(:,1)               <= this.tuCoordinates(k) & ...
						dsPosStrands(:,1) + lengths - 1 >= this.tuCoordinates(k));
					for j = 1:numel(rgnIdxs)
						thresholdedSigma = min(...
							max(sigmas(rgnIdxs(j)), this.foldChangeLowerSigmaLimit),...
							this.foldChangeUpperSigmaLimit);
					
						foldChange(...
							this.tuIndexs(k), ceil(dsPosStrands(rgnIdxs(j), 2)/2)) = ...
							this.foldChangeIntercepts(k) + this.foldChangeSlopes(k) * thresholdedSigma;
					end
				end
			elseif i == 14
				for k = 3
				rgnIdxs = find(...
						dsPosStrands(:,1)               <= this.tuCoordinates(k) & ...
						dsPosStrands(:,1) + lengths - 1 >= this.tuCoordinates(k));
					for j = 1:numel(rgnIdxs)
						thresholdedSigma = min(...
							max(sigmas(rgnIdxs(j)), this.foldChangeLowerSigmaLimit),...
							this.foldChangeUpperSigmaLimit);
					
						foldChange(...
							this.tuIndexs(k), ceil(dsPosStrands(rgnIdxs(j), 2)/2)) = ...
							this.foldChangeIntercepts(k) + this.foldChangeSlopes(k) * thresholdedSigma;
					end
				end
			end
        end
        
        function value = buildEnzymeProperties(this)
			for i = 1:16
            value = [...
                struct(...
                    'idx', this.enzymeIndexs_gyrase,...
                    'sigmaLimit', this.gyraseSigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas{i} > this.gyraseSigmaLimit, ...
                    'atpCost', this.gyraseATPCost,...
                    'deltaLK', this.gyraseDeltaLK,...
                    'activityRate', this.gyraseActivityRate,...
                    'meanDwellTime', this.gyraseMeanDwellTime,...
                    'probOfActivityFunc', @(sigma, c) 1 ./ (1 + exp(this.gyrLogisiticConst * (sigma - c.equilibriumSuperhelicalDensity))));
                struct(...
                    'idx', this.enzymeIndexs_topoIV,...
                    'sigmaLimit', this.topoIVSigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas{i} > this.topoIVSigmaLimit, ...
                    'atpCost', this.topoIVATPCost,...
                    'deltaLK', this.topoIVDeltaLK,...
                    'activityRate', this.topoIVActivityRate,...
                    'meanDwellTime', NaN,...
                    'probOfActivityFunc', @(sigma, c) ones(size(sigma)));
                struct(...
                    'idx', this.enzymeIndexs_topoI,...
                    'sigmaLimit', this.topoISigmaLimit, ...
                    'isActiveInRegion', @(sigmas) sigmas{i} < this.topoISigmaLimit, ...
                    'atpCost', this.topoIATPCost,...
                    'deltaLK', this.topoIDeltaLK,...
                    'activityRate', this.topoIActivityRate,...
                    'meanDwellTime', 0,...
                    'probOfActivityFunc', @(sigma, c) 1 ./ (1 + exp(this.topoILogisiticConst * (sigma - c.equilibriumSuperhelicalDensity))))];
		  end
	   end
    end
end

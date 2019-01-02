%CellGeometry
%
% @wholeCellModelID State_CellGeometry
% @name             Cell geometry
% @description
%
% The Cell Shape process calculates the length and surface area of the cell
% across its lifespan. It also keeps track of the cell volume. The cell is
% approximated to be of a rod shape similar to E. coli, even though
% M. genitalium is known to have a less uniform flask shape. Initially, the
% cell is modeled as a cylinder with two hemispherical caps. Once cell
% pinching commences at the midline of the cell, the shape and size of a
% "septum region" is also modeled.
%
% General geometric equations to represent the shape of a cell, and the
% idea that the density of the cell does not change across the cell cycle
% are borrowed from Domach et al. (1983). We add the assumption that the
% width of the cell remains constant during the lifespan of a cell. Thus,
% the density and cell width are inputs into our model. The cell density of
% E. coli is used 1100g/L (Baldwin et al., 1995). The cell width is
% calculated based on the initial cell mass and density and the assumption
% that the cell is a sphere. The initial cell mass is fit to result in a
% cell width of 200nm (Lind et al, 1984).
%
% Our model calculates the mass of the cell at all timesteps. The mass and
% constant density provide us with the volume at all time points.
% Volume = [2 hemispheres] + [2 cylinders] + [septum region]
% The volume of the septum region is calculated as a cylinder of length,
% 2*septumLength, and width of the cell. Then two cones (height=septum,
% radius=septum) are subtracted from this cylinder.
% (This approximation was used in Shuler et al. 1979).
%
% Since we know the volume of the cell from the cell's mass and density,
% and the septum length from our cytokinesis process, the volume formula
% gives us the length of the cell at each timestep.
%
% Similarly, we can also calculate the surface area of the cell.
% Surface Area = [2 hemispheres] + [2 cylinders] + [septum region]
% The surface area of the septum region is calculated as a cylinder of
% length, 2*septum, and width of the cell.
%
% References
% ==================
% 1. Shuler, M.L., Leung, S., Dick, C.C. (1979). A Mathematical Model for the
%    Growth of a Single Bacteria Cell. Annals of the New York Academy of Sciences
%    326: 35-52.
% 2. Domach, M.M., Leung, S.K., Cahn, R.E., Cocks, G.G., Shuler, M.L. (1983).
%    Computer model for glucose-limited growth of a single cell of Escherichia
%    coli B/r-A. Biotechnology and Bioengineering 26: 203-216.
% 3. Lind, K., Lindhardt, B., Schutten, H.J., Blom, J., Christiansen, C.
%    (1984). Serological Cross-Reactions Between Mycoplasma genitalium and
%    Mycoplasma pneumoniae. Journal of Clinical Microbiology 20: 1036-1043.
% 4. Baldwin WW, Myer R, Powell N, Anderson E, Koch AL. (1995). Buoyant
%    density of Escherichia coli is determined solely by the osmolarity of the
%    culture medium. Arch Microbiology 164: 155-157.
%
% Author: Jayodita Sanghvi, jayodita@stanford.edu
% Author: Jonathan Karr, jkarr@stanford.edu
% Author: Jared Jacobs, jmjacobs@stanford.edu
% Affilitation: Covert Lab, Department of Bioengineering, Stanford University
% Last updated: 6/3/2010
classdef CellGeometry < edu.jiangnan.fmme.cell.sim.CellState
    %property annotations
    properties (Constant)
        optionNames             = {   %names of properties that are options
            'verbosity';
            'seed';
            };
        fixedConstantNames      = {   %names of process properties that are considered fixed constants
		    'density';
			};
        fittedConstantNames     = {}; %names of process properties that are considered fitted constants, and should be stored with the simulation as such  
        stateNames              = {   %names of properties which are part of the simulation's state
            'length';
			'motherwidth';
			'daughterwidth';
			'motherlength';
			'daughterlength';
            };
        dependentStateNames     = {   %names of properties which can be calculated from the simulation's state
            'volume'
            'surfaceArea'
            %'totalLength'
            'pinched'
            'chamberVolume'
            };
    end
    
    properties (Constant)
        dryWeight = 0;      %dry weight of this class' state properties
    end

    %fixed biological constants
    properties
        density             %g/L [PUB_0553]
    end

    %state
    properties
        length               %m; width of the cell        
        motherwidth		     %m; short axis of the mother cell
		daughterwidth        %m; short axis of the daughter cell
		motherlength		 %m; long axis of the mother cell
		daughterlength		 %m; long axis of the daughter cell		
    end   
    
    %dependent state
    properties
        volume               %L
		mothervolume		 %L
		daughtervolume		 %L
    end
    
    properties (Dependent = true, SetAccess = protected)
        %cylindricalLength    %m; length of just the cylindrical part of the cell
        surfaceArea           %m^2; surface area of the cell
        pinched               %whether cytokinesis is finished
        
        chamberVolume         %volume of chamber in which cell is simulated (L)
    end
    
    %references to other parts of the cell state
    properties
        mass
    end

    %constructor
    methods
        function this = CellGeometry(wholeCellModelID, name)
            this = this@edu.jiangnan.fmme.cell.sim.CellState(wholeCellModelID, name);            
        end
    end
    
    methods
        %build object graph
        function storeObjectReferences(this, simulation)
            this.mass = simulation.state('Mass');
        end
    end
    
    %allocate memory for state
    methods
        function allocateMemory(this, numTimePoints)            
            this.motherwidth        = zeros(1, 1, numTimePoints);
            this.daughterwidth      = zeros(1, 1, numTimePoints);
            this.motherlength       = zeros(1, 1, numTimePoints);
            this.daughterlength     = zeros(1, 1, numTimePoints);			
			this.length             = zeros(1, 1, numTimePoints);
        end
    end

    %model
    methods
        %initialize to spherical shape
        function initialize(this)
            this.calculateVolume();
			[w, s] = this.calculateWidth(this.volume);
            
			this.motherwidth     = w;
			this.daughterwidth   = s;
			this.motherlength    = 1.5 * w;
			this.daughterlength  = 1.5 * s;
			this.length          = 1.5 * (w + s);
			this.mothervolume    = 4 * pi * (w/2) * (w/2) * (1.5 * w/2) / 3;
			this.daughtervolume  = 4 * pi * (s/2) * (s/2) * (1.5 * s/2) / 3;
        end
                
        %calculates width
        function [w, s] = calculateWidth(this, volume)
			mothervolume   = 3.9132e-14;
			if  mothervolume > volume
				mothervolume = volume;
				daughtervolume = 0;
			else 
				mothervolume = 3.9132e-14;
				daughtervolume = volume - mothervolume;
			end	
			w = (4 * mothervolume / pi)^(1/3);
			s = (4 * daughtervolume / pi)^(1/3);
        end
		
        %updates volume; called by metabolism
        function calculateVolume(this)
            this.volume = sum(this.mass.cell) / this.density;
        end
    end
       
    %getters
    methods
        function value = get.length(this)
            if this.daughtervolume < this.mothervolume	%0.67 * mothervolume		
                value = this.length;
            else
                value = -1;
            end
        end  
		
        function value = get.surfaceArea(this)
				 value = this.calcGeometry();
        end
		
        function result = get.pinched(this)
			result = this.daughtervolume / this.mothervolume > 2/3;
        end
		
		function surfaceArea = calcGeometry(this)
            import edu.jiangnan.fmme.cell.sim.state.CellGeometry;
            surfaceArea = ...
                CellGeometry.calculateGeometry(this.mothervolume, this.daughtervolume, this.motherwidth, this.daughterwidth);
        end
       
        function value = get.chamberVolume(this)
            value = this.mass.cellInitialDryWeight / (1 - this.mass.fractionWetWeight) / ...
                this.mass.initialBiomassConcentration;
        end 
    end
    methods (Static)
		function surfaceArea = calculateGeometry(mothervolume, daughtervolume, motherwidth, daughterwidth)
			if  daughtervolume / mothervolume == 0
				w  = motherwidth;
				s  = daughterwidth;
				%l = length;
				%surfaceArea = 2^(1/2)/2 * w * (l^2 + w^2)^(1/2);
				surfaceArea = 4 * pi * ((w/2) * (w/2) + (w/2) * 1.5 * (w/2) + (w/2) * 1.5 * (w/2));
			
			elseif daughtervolume / mothervolume < 2/3
				w = motherwidth;
				s = daughterwidth;
				%l = length;			
				%surfaceArea = 26^(1/2)/4 * w^2 + 2^(1/2)/2 * s * ((l-1.5*w)^2 + s^2)^(1/2);	
				surfaceArea = 4 * pi * ((w/2) * (w/2) + (w/2) * 1.5 * (w/2) + (w/2) * 1.5 * (w/2)) + ...
							  4 * pi * ((s/2) * (s/2) + (s/2) * 1.5 * (s/2) + (s/2) * 1.5 * (s/2));		
			else
			% cell has divided
				warning('WholeCell:warning', 'Cell has divided. Cell shape undefined');
				surfaceArea = -1;
			end
		end
    end
end

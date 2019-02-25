classdef calcRegionWeights < handle
Methods
	function weights = calcRegionWeights(lens)
	weights = max(0, lens - footprint + 1);
	end
end
end
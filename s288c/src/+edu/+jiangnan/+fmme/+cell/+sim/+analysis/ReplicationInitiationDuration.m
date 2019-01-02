%ReplicationInitiationDuration
%
% Author: Jonathan Karr, jkarr@stanford.edu
% Affiliation: Covert Lab, Department of Bioengineering, Stanford University
% Last Updated: 9/12/2010
classdef ReplicationInitiationDuration
    methods (Static)
        %inputFilePath is a path to .mat files containing replication initiation
        %durations generated by
        %edu.jiangnan.fmme.cell.sim.process.ReplicationInitiation_Test.meanDur
        %ation
        function run(inputFilePath, outputFileName)
            import edu.jiangnan.fmme.cell.sim.analysis.ReplicationInitiationDuration;
            import edu.jiangnan.fmme.cell.sim.util.PrintUtil;            
            
            %load data
            [durations, constants] = ReplicationInitiationDuration.loadDurations(inputFilePath);
            
            %statistics
            [vals, ~, means, stds, mins, maxs, meds, p10, p25, p75, p90] = ...
                ReplicationInitiationDuration.calculateStatistics(durations, constants);
            colLabels = {'Mean DnaA Copy Number', 'Site Cooperativity', 'State Cooperativity', 'No. DnaA Boxes', ...
                'Mean', 'Std. Dev', 'Min', '10%', '25%', 'Median', '75%', '90%', 'Max'};
            if nargin == 1
                PrintUtil.printToStdIO(num2cell([vals means stds mins p10 p25 meds p75 p90 maxs]), colLabels);
            else
                PrintUtil.printToFile(num2cell([vals means stds mins p10 p25 meds p75 p90 maxs]), ...
                    colLabels, [outputFileName '.xls'], 'Statistics')
            end
            
            %plots
            if nargin == 2, [axesHandle, figHandle] = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle(); end
            
            if nargin == 1, [axesHandle, figHandle] = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle(); end
            cla(axesHandle);
            tfs = ismember(constants(:, [1 3 4]), constants(1, [1 3 4]), 'rows');
            ReplicationInitiationDuration.plotScatterPlot(axesHandle, durations(tfs), constants(tfs, :));
            if nargin == 2, saveas(figHandle, [outputFileName '-ScatterPlot.pdf']); end
            
            if nargin == 1, [axesHandle, figHandle] = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle(); end
            cla(axesHandle);
            tfs = ismember(constants, constants(1, :), 'rows');
            ReplicationInitiationDuration.plotHistogram(axesHandle, durations(tfs), constants(tfs, :));
            if nargin == 2, saveas(figHandle, [outputFileName '-Histogram.pdf']); end
            
            if nargin == 1, [axesHandle, figHandle] = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle(); end
            cla(axesHandle);
            ReplicationInitiationDuration.plotDistributions(axesHandle, durations, constants);
            if nargin == 2, saveas(figHandle, [outputFileName '-Distributions.pdf']); end
            
            if nargin == 1, [axesHandle, figHandle] = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle(); end
            cla(axesHandle);
            ReplicationInitiationDuration.plotMeanVsStd(axesHandle, durations, constants);
            if nargin == 2, saveas(figHandle, [outputFileName '-MeanVsStd.pdf']); end
            
            if nargin == 2, close(figHandle); end
        end
    end
    
    methods (Static)
        function [durations, constants] = loadDurations(filePath)
            fileNames = dir([filePath filesep '*.mat']);
            durations = zeros(0, 1);
            constants = zeros(0, 4);
            path = fileparts([filePath filesep]);
            for i = 1:size(fileNames, 1)
                [~, parameters] = regexp(fileNames(i).name, ...
                    'ReplicationInitiationDuration-(?<meanDnaACopyNumber>\d+)-(?<siteCooperativity>\d+\.\d+)-(?<stateCooperativity>\d+\.\d+)-(?<nDnaABoxes>\d+)\.mat', 'match', 'names');
                
                data = load([path filesep fileNames(i).name]);
                durations = [durations; data.durations]; %#ok<AGROW>
                constants = [constants; repmat([...
                    str2double(parameters.meanDnaACopyNumber) ...
                    str2double(parameters.siteCooperativity) ...
                    str2double(parameters.stateCooperativity) ...
                    str2double(parameters.nDnaABoxes) ...
                    ], numel(data.durations), 1)]; %#ok<AGROW>
            end
            
            durations = durations / 3600;
        end
        
        function [vals, cnts, means, stds, mins, maxs, meds, p10, p25, p75, p90] = ...
                calculateStatistics(durations, constants)
            durationBins = 0:2:20;
            
            vals = unique(constants, 'rows');
            cnts = zeros(numel(durationBins), size(vals, 1));
            
            means = zeros(size(vals, 1), 1);
            stds = zeros(size(vals, 1), 1);
            mins = zeros(size(vals, 1), 1);
            maxs = zeros(size(vals, 1), 1);
            meds = zeros(size(vals, 1), 1);
            p25 = zeros(size(vals, 1), 1);
            p75 = zeros(size(vals, 1), 1);
            p90 = zeros(size(vals, 1), 1);
            p10 = zeros(size(vals, 1), 1);
            
            durations = min(durations, 3*23825/3600);
            
            for i = 1:size(vals, 1)
                tmpDurations = durations(ismember(constants, vals(i, :), 'rows'), 1);
                cnts(:,i) = histc(tmpDurations, durationBins);
                means(i) = mean(tmpDurations);
                stds(i) = std(tmpDurations);
                mins(i) = min(tmpDurations);
                maxs(i) = max(tmpDurations);
                meds(i) = median(tmpDurations);
                p25(i) = quantile(tmpDurations, 0.25);
                p75(i) = quantile(tmpDurations, 0.75);
                p90(i) = quantile(tmpDurations, 0.90);
                p10(i) = quantile(tmpDurations, 0.10);
            end
        end
    end
    
    %plots
    methods (Static)
        function plotScatterPlot(axesHandle, durations, constants)
            if isempty(axesHandle)
                axesHandle = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            plot(axesHandle, constants(:, 2), durations, '.');
            xlabel(axesHandle, 'Site Cooperativity Constant', 'fontsize', 14);
            ylabel(axesHandle, 'Replication Initiation Time (h)', 'fontsize', 14);
            title(axesHandle, sprintf('DnaA %d, State %f, DnaA Boxes %d', ...
                constants(1), constants(3), constants(4)), 'fontsize', 14);
        end
        
        function plotHistogram(axesHandle, durations, constants)
            if isempty(axesHandle)
                axesHandle = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            hist(axesHandle, durations);
            ylabel(axesHandle, 'Frequency', 'fontsize', 14);
            xlabel(axesHandle, 'Replication Initiation Time (h)', 'fontsize', 14);
            title(axesHandle, sprintf('DnaA %d, State %f, DnaA Boxes %d', ...
                constants(1), constants(3), constants(4)), 'fontsize', 14);
        end
        
        function plotDistributions(axesHandle, durations, constants)
            if isempty(axesHandle)
                axesHandle = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            [vals, cnts] = edu.jiangnan.fmme.cell.sim.analysis.ReplicationInitiationDuration.calculateStatistics(...
                durations, constants);
            
            colors = [
                1.00 0 0;
                0.75 0 0;
                0.50 0 0;
                0.25 0 0;
                0 0 0;
                0 0 0.25;
                0 0 0.50;
                0 0 0.75;
                0 0 1];
            
            h = plot(axesHandle, (0:2:20)', cnts);
            for i = 1:size(vals, 1)
                set(h(i), 'color', colors(mod(i-1, size(colors, 1))+1, :));
            end
            legend(h, cellfun(@(val1, val2, val3, val4) sprintf('%d %.2f %.2f %d', val1, val2, val3, val4), ...
                num2cell(vals(:, 1)), ...
                num2cell(vals(:, 2)), ...
                num2cell(vals(:, 3)), ...
                num2cell(vals(:, 4)), ...
                'UniformOutput', false), 'Location', 'BestOutside');
            ylabel(axesHandle, 'Frequency', 'fontsize', 14);
            xlabel(axesHandle, 'Replication Initiation Time (h)', 'fontsize', 14);
            ylim(axesHandle, [0 max(cnts(:)) + eps]);
        end
        
        function plotMeanVsStd(axesHandle, durations, constants)
            if isempty(axesHandle)
                axesHandle = edu.jiangnan.fmme.cell.sim.util.PlotUtil.newAxesHandle();
            end
            
            [vals, ~, means, stds] = edu.jiangnan.fmme.cell.sim.analysis.ReplicationInitiationDuration.calculateStatistics(...
                durations, constants);
            
            vals3 = unique(vals(:, 3));
            markerSizes = repmat([10 18 8 12 24]', ceil(numel(vals3)/5), 1);
            markerSizes = markerSizes(1:numel(vals3));
            vals2 = unique(vals(:, 2));
            markerStyles = repmat(('.+*oxsdv^<')', ceil(numel(vals2)/10), 1);
            markerStyles = markerStyles(1:numel(vals2));
            vals1 = unique(vals(:, 1));
            colors = repmat([
                1 0 0
                0 0 0
                0 0 1
                0.75 0 0
                0.50 0 0
                0.25 0 0
                0 0 0.25
                0 0 0.50
                0 0 0.75
                ], ceil(numel(vals1)/9), 1);
            colors = colors(1:numel(vals1), :);
            
            hold(axesHandle, 'on');
            for i = 1:size(vals, 1)
                h(i) = plot(axesHandle, means(i), stds(i), '.', ...
                    'MarkerSize', markerSizes(ismember(vals3, vals(i, 3))), ...
                    'Marker', markerStyles(ismember(vals2, vals(i, 2))), ...
                    'Color', colors(ismember(vals1, vals(i, 1)), :));
            end
            legend(h, cellfun(@(val1, val2, val3, val4) sprintf('%d %.2f %.2f %d', val1, val2, val3, val4), ...
                num2cell(vals(:, 1)), ...
                num2cell(vals(:, 2)), ...
                num2cell(vals(:, 3)), ...
                num2cell(vals(:, 4)), ...
                'UniformOutput', false), 'Location', 'BestOutside');
            xlim(axesHandle, [min(means) max(means)]+[-1 1]*0.1*range(means));
            ylim(axesHandle, [min(stds) max(stds)]+[-1 1]*0.1*range(stds));
            ylabel(axesHandle, 'Duration standard deviation (h)', 'fontsize', 14);
            xlabel(axesHandle, 'Mean duration (h)', 'fontsize', 14);
        end
    end
end
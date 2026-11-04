
%% Step 19a: Check sample quality with IQR and Spearman Corr
% I decided to use the RF-normalized data (chosen in Step 18).


CORR_TYPE = 'Spearman';

chosenLog  = rfLog;
iqrChosen  = iqr(chosenLog);
corrChosen = SampleCorrelation(chosenLog, CORR_TYPE);
corrChosen = corrChosen(:)';                 % turn into a row for plotting
 
scatterplotMarkOutliers(iqrChosen, 'ColumnLabels', sampleIDs, ...
    'PlotTitle', 'QC: IQR per sample (RF, log2)', ...
    'xlabel', 'Sample', 'ylabel', 'IQR');
exportgraphics(gcf, 'Fig4a_QC_IQR.png', 'Resolution', 300);
 
scatterplotMarkOutliers(corrChosen, 'ColumnLabels', sampleIDs, ...
    'PlotTitle', 'QC: average Spearman correlation (RF, log2)', ...
    'xlabel', 'Sample', 'ylabel', 'Average Spearman correlation');
exportgraphics(gcf, 'Fig4b_QC_corr.png', 'Resolution', 300);
 
 
%% Step 19b: After generating the plots, I realized that some samples have an IQR 
% of 0 and lower correlation, but they are not below the outlier cut off red line. 
% I want to double check how many miRNAs each sample really had.
 
% Percent of miRNAs with zero reads in each sample

pctZero = 100 * mean(rawCounts == 0);
disp(table(sampleIDs', pctZero', 'VariableNames', {'Sample', 'PercentZero'}))
 
% Number of miRNAs each sample had in the original ICGC file
dlbclRows = mirna_seq(ismember(mirna_seq.icgc_sample_id, sampleIDs), :);
disp(groupcounts(dlbclRows, 'icgc_sample_id'))
 
 
%% Step 19c: Remove the bad samples
% Visually, the last 6 samples (SA321182 to SA321214) look like the batch effect
% Samples only had 18-24% of the 5,850 miRNAs. The rest were empty (zero) 
% and also has an IQR of 0, fewer reads (4.7-22 million), and lower correlation (0.59-0.61, vs 0.65-0.68 for the others).

% Remove samples by their ID, not by column number.
% RF normalization treats each sample on its own, so we do not need to
% normalize again after removing samples.

REMOVE_IDS = {'SA321075', 'SA321182', 'SA321186', 'SA321190', ...
              'SA321198', 'SA321206', 'SA321214'};
 
keepSamples = ~ismember(sampleIDs, REMOVE_IDS);   % true = keep this sample
 
% Remove them from the data AND from every ID list, so everything still lines up
chosenLogQC = chosenLog(:, keepSamples);
rawCountsQC = rawCounts(:, keepSamples);
sampleIDsQC = sampleIDs(keepSamples);
donorIDsQC  = donorIDs(keepSamples);
 
fprintf('Removed %d samples; %d remain.\n', sum(~keepSamples), sum(keepSamples));


%% Step 19d: Check quality again after removing samples
% Replot the IQR amd spearman after removing 
iqrQC  = iqr(chosenLogQC);
corrQC = SampleCorrelation(chosenLogQC, CORR_TYPE);
corrQC = corrQC(:)';
 
scatterplotMarkOutliers(iqrQC, 'ColumnLabels', sampleIDsQC, ...
    'PlotTitle', 'QC after removal: IQR', ...
    'xlabel', 'Sample', 'ylabel', 'IQR');
exportgraphics(gcf, 'Fig5a_QCafter_IQR.png', 'Resolution', 300);
 
scatterplotMarkOutliers(corrQC, 'ColumnLabels', sampleIDsQC, ...
    'PlotTitle', 'QC after removal: average Spearman correlation', ...
    'xlabel', 'Sample', 'ylabel', 'Average Spearman correlation');
exportgraphics(gcf, 'Fig5b_QCafter_corr.png', 'Resolution', 300);
 
 
%% Step 20: Compare different Q values
% How MarkLowCounts works (lecture slides 134-135):
%   1. Find one cut-off value T for the whole dataset: T is the value
%      that Q of all numbers in the data fall below (e.g. Q = 0.75 -> 75%).
%   2. Keep a miRNA if it is above T in at least one sample.
%      Remove it if it is below T in every sample.
% The lecture uses Q = 0.9 by default. A lower Q (e.g. 0.75) is less strict
% and keeps more miRNAs.

Q_OPTIONS = [0.5 0.75 0.9];
nMirna = size(chosenLogQC, 1);
 
for q = Q_OPTIONS
    nLow = sum(MarkLowCounts(chosenLogQC, q));
    fprintf('Q = %.2f removes %d of %d miRNAs\n', q, nLow, size(chosenLogQC, 1));
end

%% Run the above step seperately and in the command window, 
% how much miRNAs each Q value removes should pop up
 

%% Step 20: Remove low-count miRNAs
% I chose Q = 0.75 

Q_LEVEL = 0.75;
 
isLow = MarkLowCounts(chosenLogQC, Q_LEVEL);
filteredLog      = chosenLogQC(~isLow, :);
mirnaIDsFiltered = mirnaIDs(~isLow);
fprintf('Q = %.2f: removed %d miRNAs, %d remain.\n', Q_LEVEL, sum(isLow), size(filteredLog, 1));
 

%% Step 21: Final plots of the cleaned data
% Box plot: the spread of values in each sample
figure;
boxplot(filteredLog, 'Labels', sampleIDsQC, 'LabelOrientation', 'inline');
title('Final data: RF-normalized, log2, low counts removed');
ylabel('log2(RF-normalized count)');
exportgraphics(gcf, 'Fig6a_final_boxplot.png', 'Resolution', 300);
 
% IQR per sample
iqrFinal = iqr(filteredLog);
scatterplotMarkOutliers(iqrFinal, 'ColumnLabels', sampleIDsQC, ...
    'PlotTitle', 'Final data: IQR per sample', ...
    'xlabel', 'Sample', 'ylabel', 'IQR');
exportgraphics(gcf, 'Fig6b_final_IQR.png', 'Resolution', 300);
 
% Average correlation per sample
corrFinal = SampleCorrelation(filteredLog, CORR_TYPE);
corrFinal = corrFinal(:)';
scatterplotMarkOutliers(corrFinal, 'ColumnLabels', sampleIDsQC, ...
    'PlotTitle', 'Final data: average Spearman correlation', ...
    'xlabel', 'Sample', 'ylabel', 'Average Spearman correlation');
exportgraphics(gcf, 'Fig6c_final_corr.png', 'Resolution', 300);
 
 

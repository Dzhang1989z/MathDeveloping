clear all; clc; close all;

%% sub part 
clear all;
InputGeneTotal = 282;
CategoryGeneTotal = 20720;

[~, ~, alldata] = xlsread('GOanalysis_SUB_NORM_NEW.xlsx');
GeneRatio = cell2mat(alldata(2:end,7))/InputGeneTotal;
BgRatio = cell2mat(alldata(2:end,8))/CategoryGeneTotal;

FoldEnrichment = log2(GeneRatio./BgRatio);
Log10Pvalue = -log10(cell2mat(alldata(2:end,3)));
DotSize = cell2mat(alldata(2:end,7));

figure(2);
Idx = min(find(cell2mat(alldata(2:end,6))>0.05));
for p = 1:(Idx-1)
    scatter(FoldEnrichment(p), Log10Pvalue(p), DotSize(p),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor','c','MarkerFaceAlpha',1);
    hold on;
end
for p = Idx:length(FoldEnrichment)
    scatter(FoldEnrichment(p), Log10Pvalue(p), DotSize(p),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor','c','MarkerFaceAlpha',0.2);
    hold on;
end
% xlim([1, 3]);
ylim([1, 6]);
xlabel('Fold enrichment ratios (log2)');
ylabel('P value (-log10)');
set(gca, 'FontSize', 12);
set(gca,'LineWidth',1.5);


%% mul part 
clear all;

InputGeneTotal = 660;
CategoryGeneTotal = 20720;

[~, ~, alldata] = xlsread('GOanalysis_mul_NORM_NEW.xlsx');
GeneRatio = cell2mat(alldata(2:end,7))/InputGeneTotal;
BgRatio = cell2mat(alldata(2:end,8))/CategoryGeneTotal;

FoldEnrichment = log2(GeneRatio./BgRatio);
Log10Pvalue = -log10(cell2mat(alldata(2:end,3)));
DotSize = cell2mat(alldata(2:end,7));

figure(3);
Idx = min(find(cell2mat(alldata(2:end,6))>0.05));
for p = 1:(Idx-1)
    scatter(FoldEnrichment(p), Log10Pvalue(p), DotSize(p),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor','r','MarkerFaceAlpha',1);
    hold on;
end
for p = Idx:length(FoldEnrichment)
    scatter(FoldEnrichment(p), Log10Pvalue(p), DotSize(p),'MarkerEdgeColor','k','MarkerEdgeAlpha',0,'MarkerFaceColor','r','MarkerFaceAlpha',0.2);
    hold on;
end
% xlim([1, 3]);
ylim([1, 15]);
xlabel('Fold enrichment ratios (log2)');
ylabel('P value (-log10)');
set(gca, 'FontSize', 12);
set(gca,'LineWidth',1.5);

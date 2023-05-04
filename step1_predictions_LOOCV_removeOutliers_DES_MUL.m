clear all; clc; close all;

%% build T1 prediction model
% volume data includes 3 measures
% including SurfArea, GrayVol, ThickAvg is 132x68 each; 
load('DES_T1_VolumeData.mat');  
% behavior data includes 9 performance.
% all 9 behaviors is 132x1
%  behavior data including: 'CMAT_Addition_StS', 'CMAT_Subtraction_StS', 'CMAT_Multiplication_StS', 'CMAT_Division_StS', 'CMAT_BasicCalc_Comp_Quotient', 'CMAT_SUM','age_ses-T1', 'gender', 'eTIV'
load('T1_BehaviorData.mat'); %% 
behav_vec = BehaviorData.CMAT_Multiplication_StS;
behav_vec2 = BehaviorData.CMAT_Subtraction_StS;

% preparing data
GrayVol = VolumeData.GrayVol;     %  132x68 

AgeT1 = BehaviorData.AGE;
eTIVT1 = BehaviorData.eTIV;
IQT1 =  cell2mat(BehaviorData.IQ_FullScale);

%% 去除GMV异常的被试
dropIdxArray = [26 105];    % remove two outliers

GrayVol(dropIdxArray,:) = [];
behav_vec(dropIdxArray) = [];
behav_vec2(dropIdxArray) = [];

AgeT1(dropIdxArray) = [];
eTIVT1(dropIdxArray) = [];
IQT1(dropIdxArray) = [];

%% LOOCV 查看GMV异常被试
numROIs = size(GrayVol, 2);
thresh = 0.01;
WeightArray = zeros(size(GrayVol,1), size(GrayVol,2)+3);
for leftout = 1:size(GrayVol,1)

    trainIdx = setdiff([1:size(GrayVol,1)],leftout);
    TrainingData = GrayVol(trainIdx,:);
    TrainY = behav_vec(trainIdx)';
    TestingData = GrayVol(leftout,:);
    TestingY = behav_vec(leftout)';
    
    r_mat = []; p_mat = [];    
    for roiIdx = 1:size(GrayVol,2)
        [r, p] = partialcorr(TrainingData(:,roiIdx), TrainY, [AgeT1(trainIdx)', eTIVT1(trainIdx)', behav_vec2(trainIdx)' ]);
        r_mat(roiIdx) = r;
        p_mat(roiIdx) = p;
    end
    pos_mask = zeros(numROIs, 1);
    neg_mask = zeros(numROIs, 1);
    pos_nodes = find(p_mat < thresh);

    % train linear regression model based on the volume feature, just positive
    % features
    TrainX = [TrainingData(:,pos_nodes), AgeT1(trainIdx)', eTIVT1(trainIdx)', ones(length((trainIdx)),1)];  % 132x11
    [Beta, BINT, R, RINT, stats] = regress(TrainY, TrainX);
    
    % save weight 
    WeightArray(leftout, pos_nodes) = Beta(1:end-3);
    if length(Beta) > 3
        WeightArray(leftout, end-3:end) = Beta(end-3:end);
    end

    % testing part 
    TestX = [TestingData(:,pos_nodes), AgeT1(leftout)', eTIVT1(leftout)',  ones(1,1)];  % 132x11
    TestingY_pred = TestX*Beta;

    behav_vec_pred(leftout) = TestingY_pred;
end

MUL_GMV_prediction = [behav_vec', behav_vec_pred'];    % SUB_GMV_prediction = 130x2
csvwrite('MUL_GMV_prediction.csv',MUL_GMV_prediction);

% plot parts
figure(1);
[R_pos, P_pos] = corr(behav_vec', behav_vec_pred');
x = behav_vec;
y = behav_vec_pred;
[h, sortorder] = sort(x);
x = x(sortorder);
y = y(sortorder);
t = figure(1);
plot(x, y, 'bx');
hold on;
[p, s] = polyfit(x, y, 1);
[yfit, dy] = polyconf(p, x, s, 'predopt', 'curve');
plot(x, yfit, 'color', 'r');
plot(x, yfit-dy, 'color', 'r','linestyle',':');
plot(x, yfit+dy, 'color', 'r','linestyle',':');
xlabel('Observed preformance');
ylabel('Predictied preformance');
title({'Predicition model', ['r = ' num2str(R_pos) ', p = ' num2str(P_pos) ' on testing samples']}, 'FontSize', 12);

% permutation test
R_pos_array = [];
for p = 1:10000
    permIdx = randperm(length(behav_vec));
    behav_vec_Perm = behav_vec(permIdx);
    [R_pos_Perm, drop] = corr(behav_vec_pred', behav_vec_Perm');
    R_pos_array(p) = R_pos_Perm;
end
% 
R_pos_model = R_pos;
permP = length(find(R_pos_array>R_pos_model))/length(R_pos_array);
figure(2);
hist(R_pos_array, 24);
hold on;
plot([R_pos_model R_pos_model], [0 1400], 'r-');
xlabel('R value');
ylabel('Frequency');
title(['Permutation test p = ' num2str(permP)])

%% testing on T2 samples
WeightArray = mean(WeightArray,1);
save('WeightArray_mul_removeOutliers','WeightArray');


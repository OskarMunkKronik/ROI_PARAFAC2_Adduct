
cd([BaseDirectory])
load([BaseDirectory,'\DataFiles\MatFiles\241101_PARAFAC2_Models.mat'])
% load([BaseDirectory,'\DataFiles\MatFiles\241010_mzroi_aug.mat'])
% load([BaseDirectory,'\DataFiles\MatFiles\241010_MSroi_aug.mat'])

%% Inputs
Options_AdductFinder.nMax = 10;
% mzroi = mzroi_aug_Int;
mzHydrogen = 1.007825031900000;
Options_AdductFinder.Adducts = {'Na','K','NH4','K-Na','Na-NH4','K-NH4'};
Options_AdductFinder.Adducts_Exact_mz = [22.98976928,38.96370648,18.034374132]-mzHydrogen; %Na, K, NH4, 
Options_AdductFinder.Adducts_Exact_mz(4) = diff(Options_AdductFinder.Adducts_Exact_mz([1,2]));
Options_AdductFinder.Adducts_Exact_mz(5) = diff(Options_AdductFinder.Adducts_Exact_mz([3,1]));
Options_AdductFinder.Adducts_Exact_mz(6) = diff(Options_AdductFinder.Adducts_Exact_mz([3,2]));



Options_AdductFinder.mzDev  = 0.01;
Options_AdductFinder.CorrThresh = 0.9
Options_AdductFinder.nMax = 3
Options_AdductFinder.RtDev = 0.05;

% %%
% %%
% col = 'rbgckym'
% F   = 12;
% Int = 42;
% 
% clf 
% load(['DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat'])
% mzroi = mzroi_aug_Int;
% Rt = Rt_Int;
% h = 1
% fVec = [1:F]
% for f = fVec
%     subplot(F,2,h)
%     for k = 1:size(Model{Int,F}.Loads{2},3)
%     plot(Model{Int,F}.Loads{2}(:,f,k).*Model{Int,F}.Loads{3}(k,f))%,'Color',col(f))
%     hold on
%     title(num2str(f))
%     end
%     h = h+2;
% end
% h = 2;
% for f = fVec
%     subplot(F,2,h)
% stem(mzroi,Model{Int,F}.Loads{1}(:,f))%,'Color',col(f))
% hold on
% h = h +2;
% end 
% axis tight
% %%
% mzTarget = [421.209, 399.2270557519000];%608.35
% clf
% for nMZ = 1:length(mzTarget)
%  [~, a] = min(abs(mzroi_aug-mzTarget(nMZ)))
% for k = 1:size(MSroi_aug,1)
%     plot(MSroi_aug{k}(a,:),'Color',col(nMZ))
%     hold on 
% end 
% end 
%%  Find adducts and flag them 
[mzInd,out,groupAdduct] = deal(cell(size(Model)));
Adduct_flag = zeros(size(Model));
for Int = 1:size(Model,1)
fprintf(1,'Interval: %i\n',Int)
load(['DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat'])
for f = 1:size(Model,2)

[mzInd{Int,f},out{Int,f},groupAdduct{Int,f}] = FindAdducts(Model{Int,f},X,mzroi_aug_Int,Rt_Int,Options_AdductFinder);
Adduct_flag(Int,f) = nnz(groupAdduct{Int,f});
end 
end 
save([BaseDirectory,'\DataFiles\MatFiles\AdductFlagging.mat'],'Adduct_flag','mzInd','out','groupAdduct','Options_AdductFinder')
% %% Recalculate PF2 models on summed data 
% [mzroi_new,X_new] = SumAdducts(mzroi_aug_Int,X,groupAdduct{Int,f});

%%

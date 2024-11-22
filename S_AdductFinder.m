
cd([BaseDirectory])
load([BaseDirectory,'\DataFiles\MatFiles\241101_PARAFAC2_Models.mat'])
%% Inputs
mzHydrogen = 1.007825031900000;
Options_AdductFinder.Adducts = {'Na','K','NH4','K-Na','Na-NH4','K-NH4'};
Options_AdductFinder.Adducts_Exact_mz = [22.98976928,38.96370648,18.034374132]-mzHydrogen;     % Na, K, NH4,
Options_AdductFinder.Adducts_Exact_mz(4) = diff(Options_AdductFinder.Adducts_Exact_mz([1,2])); % diff [M+K]  - [M+Na],
Options_AdductFinder.Adducts_Exact_mz(5) = diff(Options_AdductFinder.Adducts_Exact_mz([3,1])); % diff [M+Na] - [M+NH4+],
Options_AdductFinder.Adducts_Exact_mz(6) = diff(Options_AdductFinder.Adducts_Exact_mz([3,2])); % diff [M+K]  - [M+NH4+],
Options_AdductFinder.mzDev       = 0.01;        %allowed m/z deviation
Options_AdductFinder.CorrThresh  = 0.9;         %Elution profile correlation
Options_AdductFinder.nMax        = 3;           %number of m/z peaks to consider
Options_AdductFinder.RtDev       = 0.05;        %minutes +/-

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

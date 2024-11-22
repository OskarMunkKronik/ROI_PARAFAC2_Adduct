load([BaseDirectory,'\DataFiles\MatFiles\AdductFlagging.mat'])
load([BaseDirectory,'\DataFiles\MatFiles\Interval.mat'])

%% Make vectors ready for parrallel computing


% Interval_Vec = repelem(1:size(Intervals,1),Options.MaxFactor);
% Factor_Vec = repmat(1:Options.MaxFactor,1,size(Intervals,1));
[Interval_Vec,Factor_Vec] = find(Adduct_flag);
[Model,ProcessTime] = deal(cell(size(Intervals,1),Options.MaxFactor));
[Model_tmp,ProcessTime_tmp] = deal(cell(length(Interval_Vec),1));
%% Save new adduct preprocessed
uInt =  unique(Interval_Vec);
for nInt =1:length(uInt)
    Int = uInt(nInt);
    fprintf(1,'Interval: %i\n',Int)
    
Adduct_cat = [];
for f = 1:size(Adduct_flag,2)
    Adduct_cat = cat(1,Adduct_cat,groupAdduct{Int,f});
end 
Adduct_cat = unique(Adduct_cat,'rows');


Data = load([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat']);
[Data.mzroi_aug_Int,Data.X] = SumAdducts(Data.mzroi_aug_Int,Data.X,Adduct_cat,[mzHydrogen,mzHydrogen*2]);
mzroi_aug_Int = Data.mzroi_aug_Int;
Rt_Int        = Data.Rt_Int;
X             = Data.X;

save([BaseDirectory,'DataFiles\MatFiles\Intervals_SummedAdduct\Interval',num2str(Int),'.mat'],'mzroi_aug_Int','Rt_Int','X','Adduct_cat')
end 
%%
% waitbar
NbrePts = length(Interval_Vec);
close all force
% Waitbar's Message Field
Msg = ['PARAFAC2 Progress...!'];
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);

parfor n = 1:length(Interval_Vec)
    hWaitbarMsgQueue.send(0);
    Data = load([BaseDirectory,'DataFiles\MatFiles\Intervals_SummedAdduct\Interval',num2str(Interval_Vec(n)),'.mat']);
    [Model_tmp{n},ProcessTime_tmp{n}] = pf2flex(Data.X,Factor_Vec(n),struct('MaxIter',1500,'ConvCrit',1e-8));
end

for  n = 1:length(Interval_Vec)
    Model{Interval_Vec(n),Factor_Vec(n)} = Model_tmp{n};
    ProcessTime{Interval_Vec(n),Factor_Vec(n)} = ProcessTime_tmp{n};

end
%Commented to not overwrite the file 
% save([BaseDirectory,'DataFiles\MatFiles\','241101_PARAFAC2_Models_AdductRecalc.mat'],'Model','ProcessTime')

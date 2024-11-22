%%Import data and sort fileList and sampleList
%Sample Information
cd('I:\SCIENCE-PLEN-ECP-AnalytChem\People\Oskar\MasterThesis2\MatLab\Workflow\Import_Data')
InputData
Addpaths
cd('I:\SCIENCE-PLEN-ECP-AnalytChem\People\Oskar\MasterThesis2\DataFiles\CSV')
S_Import_SampleInformation
% % Loading MS trace
cd([BaseDirectory])

cd([BaseDirectory,'\DataFiles\NetCDF'])
fileList01  = dir('*01.CDF');
fileList02  = dir('*02.CDF');
[fileList01]=SortFileList(fileList01);
[fileList02]=SortFileList(fileList02);
%
fileList01=fileList01(SampleList{:,4}>0);
fileList02=fileList02(SampleList{:,4}>0);

nSamples=length(fileList01);
%Initialize PARAFAC2
% factor=16;
% nSpl=3;
% nGen=1;
% All=0;
% SaveFile=1;
% PLS=0;
%Sorts the sampleList and file List so that samples from the same WWTP are
%together and QC samples etc.

[SampleList2,fileList01,fileList02,SampleType] = SortSampleList2(SampleList2,fileList01,fileList02,nSamples,3);

cd D:\data\Oskar\ShortCommunication\CodeGit\roi-paper
ProjectPBA
addpath D:\data\Oskar\Pulsed2DLC\ChromTools\chromCDF
%Initialize ROI parameters
Vecmzerror= [8:4:60];%[8*5*10^-6*500];
ppm=1;
nScan=2020;
minroi=8;
thresh = 200;
sScan=1;
eScan=nScan;

% Options
Options.NumTrace = 2;
Options.MaxFactor = 12;

% Pre-allocation
[mzroi,MSroi,Rt] = deal(cell(length(fileList02),Options.NumTrace));

BaseDirectory = 'D:\data\Oskar\AdductCurveResolutionPaper\';
cd(BaseDirectory)
ProjectPBA


%% Optimize mzerror
clf
cd('I:\SCIENCE-PLEN-ECP-AnalytChem\People\Oskar\MasterThesis2\DataFiles\NetCDF\')
MassAccuracy = 10;
MZmultiplyFactors = 1:0.3:4;
% [mzroi, MSroi, RawData,runTime,sizeMZRoI,Rt,trace] = deal(cell(2,1));
OptionsROI = struct('minroi',minroi,...
    'mzerror',8,...
    'ppm',true,...
    'thresh',200,...
    'wmean',true,...
    'GapAllowed',1,...
    'NumTrace',2,...
    'verbose',true,...
    'prefilter',true,...
    'fillIn',0, ...
    'IMS',false);
save([BaseDirectory,'DataFiles\MatFiles\','241010_InputData'],'Options',"OptionsROI")
kVec  = 1:20:length(fileList01);
[mzerror_Sample,nmz] = deal(zeros(length(MZmultiplyFactors),length(kVec)));
for n = 1:length(kVec)
    k = kVec(n);
    [mzerror_Sample(:,n),~,nmz(:,n)]=OptParamRoi(fileList01(k).name,MassAccuracy,MZmultiplyFactors,OptionsROI.minroi,OptionsROI);
end
%%
cd('I:\SCIENCE-PLEN-ECP-AnalytChem\People\Oskar\MasterThesis2\DataFiles\NetCDF\')

NbrePts = length(fileList01)*2;
mzerror = median(mzerror_Sample,'all')
close all
% Waitbar's Message Field
Msg = ['Sample ',num2str(1),' - ROI Progress...!'];
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);
for T = 1:Options.NumTrace
    fprintf(1,'Trace: %i/%i\n',T,Options.NumTrace)

    parfor k = 1:length(fileList01)
        if T == 1
            FileName = fileList01(k).name;
        else 

            FileName = fileList02(k).name;
        end 
        hWaitbarMsgQueue.send(0);
              
        % ROI
        [mzroi{k,T},MSroi{k,T},~,     ~,      ~,     Rt,~,~]=ROIpeaks_ACN(FileName,struct('thresh',thresh,'mzerror',mzerror,'minroi',minroi,'ppm',ppm,'wmean',1,'prefilter',true,'GapAllowed',1,'verbose',false,'CollapseROIs',true));%,thresh,mzerror,minroi,sScan,eScan,ppm,wmean)
        % [mzroi,     MSroi,   PeaksMat,runTime,sizeMZRoI,Rt,PeaksMat1,Dt]
    end

end
OptionsROI = struct('thresh',thresh,'mzerror',mzerror,'minroi',minroi,'ppm',ppm,'wmean',1,'prefilter',true,'GapAllowed',1,'verbose',false,'CollapseROIs',true);

   [~,~,~,     ~,      ~,     Rt,~,~]=ROIpeaks_ACN(fileList01(k).name,struct('thresh',thresh,'mzerror',mzerror,'minroi',minroi,'ppm',ppm,'wmean',1,'prefilter',true,'GapAllowed',1,'verbose',false,'CollapseROIs',true));%,thresh,mzerror,minroi,sScan,eScan,ppm,wmean)
      
save([BaseDirectory,'DataFiles\MatFiles\','241010_mzroi'],'mzroi',"OptionsROI")
save([BaseDirectory,'DataFiles\MatFiles\','241010_MSroi.mat'],'MSroi',"Rt","OptionsROI")
%% Augment
runTime_aug = tic;
OptionsROI.mzerror = mzerror;
% Options = struct('thresh',thresh,'mzerror',Vecmzerror,'minroi',minroi,'ppm',false,'wmean',1,'prefilter',true,'GapAllowed',GapAllowed,'verbose',false,'CollapseROIs',true);
[mzroi_aug,MSroi_aug]=MSroiaug_ACN2(mzroi,MSroi,OptionsROI);
runTime_aug =toc(runTime_aug );

save([BaseDirectory,'DataFiles\MatFiles\','241010_mzroi_aug'],'mzroi_aug')
save([BaseDirectory,'DataFiles\MatFiles\','241010_MSroi_aug.mat'],'MSroi_aug','Rt')

%% Interval selection

%Makes mean BPC of all samples
x = zeros(nSamples,size(MSroi_aug{1},2));
for k = 1:size(MSroi_aug,1)
    x(k,:) = max(MSroi_aug{k}(:,1:2019));
end
x = mean(x,1);

%Calculates intervals
% x: The vector on which the savitsky-golay needs to be calculated
% size(1 x Number of Scans
% m: number of filtered points (FIR order)
% M: output point , M = (m+1)/2; % middle point for odd m
% n: approximation polynomial order
Intervals = Int_calc(x,3,2,12);
save([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval.mat'],'Intervals')
plot(x)
xline(Intervals(:,2))
%% saves intervals (only nnz m/zs)

for Int = 1:size(Intervals,1)
    % X = cell(size(MSroi_aug));
    X = Sparse2tensor(MSroi_aug,Intervals(Int,:));
    Ind_mz = (sum(X,[2,3])>0);
    X = X(Ind_mz,:,:);
    mzroi_aug_Int = mzroi_aug(Ind_mz);
    Rt_Int = Rt(Intervals(Int,1):Intervals(Int,2))./60;
    save([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat'],'X',"mzroi_aug_Int","Ind_mz","Rt_Int")
end

%% PARAFAC2
load('D:\data\Oskar\AdductCurveResolutionPaper\DataFiles\MatFiles\Intervals\Interval.mat')
Interval_Vec = repelem(1:size(Intervals,1),Options.MaxFactor);
Factor_Vec = repmat(1:Options.MaxFactor,1,size(Intervals,1));
[Model,ProcessTime] = deal(cell(size(Intervals,1),Options.MaxFactor));
[Model_tmp,ProcessTime_tmp] = deal(cell(length(Interval_Vec),1));

% waitbar
NbrePts = length(Interval_Vec);
close all force
% Waitbar's Message Field
Msg = ['PARAFAC2 Progress...!'];
% Create ParFor Waitbar
[hWaitbar,hWaitbarMsgQueue]= ParForWaitbarCreateMH(Msg,NbrePts);

parfor n = 1:length(Interval_Vec)
    hWaitbarMsgQueue.send(0);
    Data = load([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Interval_Vec(n)),'.mat']);
    % [Model,ProcessTime] = pf2flex(X,Factor_Vec(n),struct('MaxIter',2,'ConvCrit',1e-8));
    [Model_tmp{n},ProcessTime_tmp{n}] = pf2flex(Data.X,Factor_Vec(n),struct('MaxIter',1500,'ConvCrit',1e-8))
end

for  n = 1:length(Interval_Vec)
    Model{Interval_Vec(n),Factor_Vec(n)} = Model_tmp{n};
    ProcessTime{Interval_Vec(n),Factor_Vec(n)} = ProcessTime_tmp{n};

end
save([BaseDirectory,'DataFiles\MatFiles\','241101_PARAFAC2_Models.mat'],'Model','ProcessTime')
% clear MSroi_aug MSroi mzroi_aug mzroi
%% Plot
%Plots the elution profile
Int = Interval_Vec(100);
F = Factor_Vec(100);
Model_tmp= Model{Int,F};
load([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat']);
h = 1;
for f = 1:F
    subplot(F,3,h)
    for k = 1:size(Model_tmp.Loads{2},3)
        
        y = Model_tmp.Loads{2}(:,f,k).*Model_tmp.Loads{3}(k,f);
        plot(Rt_Int,y,'b')%,'color',color{k})
        hold on
    end
    axis tight
    xlabel('Rt (min)')
    ylabel('Intensity')

    h = h+1;
    subplot(F,3,h)
    stem(mzroi_aug_Int,Model_tmp.Loads{1}(:,f),'Marker','none')
    xlabel('m/z')
    [BasePeak,mzBasePeak]=sortrows(Model_tmp.Loads{1}(:,f),'descend');
    text(mzroi_aug_Int(mzBasePeak(1:5)),double(BasePeak(1:5)),num2str(round(mzroi_aug_Int(mzBasePeak(1:5)),3)))

    h = h+1;
    subplot(F,3,h)
    stem(Model_tmp.Loads{3}(:,f),'Marker','none')
    xlabel('sample number')
    h = h+1;

end

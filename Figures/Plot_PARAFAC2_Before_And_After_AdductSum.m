
Int = 42;
F = 5;

for n = 1:length(Interval_Vec)

    close all
    F = Factor_Vec(n)
    Int = Interval_Vec(n)

    if Int > 20 && Int < 60
        load([BaseDirectory,'\DataFiles\MatFiles\241101_PARAFAC2_Models_AdductRecalc.mat'])
        Data=load([BaseDirectory,'DataFiles\MatFiles\Intervals_SummedAdduct\Interval',num2str(Int),'.mat']);
        % Data.mzroi_aug_Int(groupAdduct{Int,F}(:,1))
        %%
        % [Data.mzroi_aug_Int,Data.X] = SumAdducts(Data.mzroi_aug_Int,Data.X,groupAdduct{Int,F},[mzHydrogen,mzHydrogen*2]);
        figure
        sgtitle(['After Summed Adduct, Interval No: ',num2str(Int)])
        h=1;
        for f = 1:F
            subplot(F,2,h)
            for k = 1:size(Model{Int,F}.Loads{2},3)
                plot(Data.Rt_Int,    Model{Int,F}.Loads{2}(:,f,k).*Model{Int,F}.Loads{3}(k,f))%,'Color',col(f))
                hold on
                title(num2str(f))
            end
            h = h+2;
        end
        h = 2;
        for f = 1:F
            subplot(F,2,h)
            stem(Data.mzroi_aug_Int,Model{Int,F}.Loads{1}(:,f),'marker','none')%,'Color',col(f))
            hold on
            h = h +2;
            xlim([50 1000])
        end
        axis tight

        %%
        load([BaseDirectory,'\DataFiles\MatFiles\241101_PARAFAC2_Models.mat'])
        h=1;
        % Int = 42;
        % F = 5;

        Data= load([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat']);
        % Data.mzroi_aug_Int(groupAdduct{Int,F}(:,1))
        %%
        % [Data.mzroi_aug_Int,Data.X] = SumAdducts(Data.mzroi_aug_Int,Data.X,groupAdduct{Int,F},[mzHydrogen,mzHydrogen*2]);
        figure
        sgtitle(['Before Summed Adduct, Interval No: ',num2str(Int)])
        h=1;
        for f = 1:F
            subplot(F,2,h)
            for k = 1:size(Model{Int,F}.Loads{2},3)
                plot(Data.Rt_Int,  Model{Int,F}.Loads{2}(:,f,k).*Model{Int,F}.Loads{3}(k,f))%,'Color',col(f))
                hold on
                title(num2str(f))
            end
            h = h+2;
        end
        h = 2;
        for f = 1:F
            subplot(F,2,h)
            stem(Data.mzroi_aug_Int,Model{Int,F}.Loads{1}(:,f),'marker','none')%,'Color',col(f))
            hold on
            h = h +2;
            xlim([50 1000])
        end
        axis tight
        pause
    end
end
% 
% function [p] = PlotPF2(mzroi,Model)
% h=1;
%         for f = 1:F
%             subplot(F,2,h)
%             for k = 1:size(Model{Int,F}.Loads{2},3)
%                 plot(Data.Rt_Int,    Model{Int,F}.Loads{2}(:,f,k).*Model{Int,F}.Loads{3}(k,f))%,'Color',col(f))
%                 hold on
%                 title(num2str(f))
%             end
%             h = h+2;
%         end
%         h = 2;
%         for f = 1:F
%             subplot(F,2,h)
%             stem(Data.mzroi_aug_Int,Model{Int,F}.Loads{1}(:,f),'marker','none')%,'Color',col(f))
%             hold on
%             h = h +2;
%             xlim([50 1000])
%         end
%         axis tight
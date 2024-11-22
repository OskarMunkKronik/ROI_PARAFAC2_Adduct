clf
%Before pre_processing
Int = 42;
Data = load([BaseDirectory,'DataFiles\MatFiles\Intervals\Interval',num2str(Int),'.mat']);
stem(Data.mzroi_aug_Int,sum(Data.X,[2,3]),'Marker','none','Color','b')
%Ater pre_processing
hold on 
load([BaseDirectory,'DataFiles\MatFiles\Intervals_SummedAdduct\Interval',num2str(Int),'.mat']);
stem(Data.mzroi_aug_Int,sum(Data.X,[2,3]),'Marker','none','Color','r')
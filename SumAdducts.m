function [mzroi,X] = SumAdducts(mzroi,X,groupAdduct,Isotopes)
% [~,~,uG] = unique(groupAdduct(:,2));
mzInd = groupAdduct(:,1);
for n = 1:max(groupAdduct(:,2))
    uG           = find(groupAdduct(:,2) == n);
    % a = [];
    for nIsotope = 1:length(Isotopes)
    [~,a] = min(abs(mzroi-mzroi(groupAdduct(uG,1))'-Isotopes(nIsotope)));
    mzInd = cat(2,mzInd,a');
    end 
    for nIsotope = 1:size(mzInd,2)
    X(mzInd(uG(1),nIsotope),:,:) = sum(X(mzInd(uG,nIsotope),:,:),1);
    X(mzInd(uG(2:end),nIsotope),:,:) = 0;
    mzroi(mzInd(uG(2:end),nIsotope)) = 0;
    end
end
Ind = sum(X,[2,3])==0;
X(Ind,:,:) = [];
mzroi(Ind) = [];
end
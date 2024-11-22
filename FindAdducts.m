
function [mzInd,out_Cell,groupAdduct] = FindAdducts(Model,X,mzroi,Rt,Options)
out_Cell = cell(length(Options.Adducts),1);
mzInd = [];
groupAdduct = [];
F = size(Model.Loads{1},2);
for nAdduct = 1:length(Options.Adducts)
    mzAdduct = Options.Adducts_Exact_mz(nAdduct);
    out = struct;
    %% Check elution profiles
    out.EluCorr = [];
    D_dist = tril(corr(squeeze(sum(Model.Loads{2}(:,:,:),3))),-1) > Options.CorrThresh;
    [Elu1,Elu2] = find(D_dist);
    if ~isempty(Elu2)
        out.EluCorr = [Elu1,Elu2];
    end
    %% Diff m/z
    out.mzAdduct = [];
    out.RtMax    = [];
    % function [IndMZ]
    % Options.nMax
    if ~isempty(out.EluCorr)
        IndMZ = zeros(size(Model.Loads{1},1),F);
        for f = 1:F
            [~,IndMZ(:,f)] = sort(Model.Loads{1}(:,f),'descend');
        end
        IndMZ = IndMZ(1:Options.nMax,:);
        for f_count = 1:length(out.EluCorr(:,1))
            f = out.EluCorr(f_count,1);
            for ff_count = 1:length(out.EluCorr(:,2))
                ff = out.EluCorr(ff_count,2);
                D_dist = abs(tril(mzroi(IndMZ(:,f)) - mzroi(IndMZ(:,ff))',-1));
                [r,c] = find(D_dist >  mzAdduct-Options.mzDev & D_dist <  mzAdduct+Options.mzDev);
                mzF1 = IndMZ(r,f);
                mzF2 = IndMZ(c,ff);
                if ~isempty(out.mzAdduct)
                    out.mzAdduct = cat(1,out.mzAdduct,[mzF1,mzF2,repelem(f,length(mzF1))',repelem(ff,length(mzF2))']);
                    % [~, RtMax(1) ]   = max(squeeze(sum(Model.Loads{2}(:,f,:),3)))
                    % [~, RtMax(2) ]   = max(squeeze(sum(Model.Loads{2}(:,ff,:),3)))
                    % out.RtMax = cat(1,out.RtMax,RtMax);

                else
                    out.mzAdduct = [mzF1,mzF2,repelem(f,length(mzF1))',repelem(ff,length(mzF2))'];
                    % [~, out.RtMax(:,1) ]   = max(squeeze(sum(Model.Loads{2}(:,f,:),3)))
                    % [~, out.RtMax(:,2) ]   = max(squeeze(sum(Model.Loads{2}(:,ff,:),3)))
                end
            end
        end

    end
    RtMax = zeros(size(out.mzAdduct,1),2);
    for n = 1:size(RtMax,1)
        [~,a]              = max(X(out.mzAdduct(n,1:2),:,:),[],2);
        RtMax(n,1:2)       = round(median(a,3));
        median(a,3)
        % [~, RtMax(n,1) ]   = max(squeeze(sum(Model.Loads{2}(:,f,:),3)));
        % [~, RtMax(n,2) ]   = max(squeeze(sum(Model.Loads{2}(:,ff,:),3)));
    end

    % RtMax = abs(diff(Rt(RtMax),[],2))  ;
    out.mzAdduct = cat(2,out.mzAdduct,RtMax);
    out_Cell{nAdduct} = out;
    if ~isempty(RtMax)
        if isempty(mzInd)
            mzInd = out.mzAdduct(:,[1,2,5,6]);
        else
            mzInd = cat(1,mzInd,out.mzAdduct(:,[1,2,5,6]));
        end


        mzInd = unique(mzInd,'rows');
        RtMax_2 = zeros(size(mzInd,1),2);
        RtMax_2(:,1:2) = Rt(mzInd(:,3:4));
        Ind = abs(diff(RtMax_2,[],2))< Options.RtDev;
        mzInd   = mzInd(Ind,:);
        RtMax_2 = RtMax_2(Ind,:);
        [~,~,uG] = unique(abs(diff(RtMax_2,[],2))< Options.RtDev);
        u = zeros(length(unique(mzInd(uG==1,1:2))),1);
        if ~isempty(mzInd)
            u(:,1) = unique(mzInd(uG==1,1:2));
            groupAdduct = cat(2,u,repelem(uG(1),length(unique(mzInd(uG==1,1:2))))');
            if max(uG)>1
                for n = 2:length(uG)
                    groupAdduct = cat(1,groupAdduct,cat(2,unique(mzInd(uG==n,1:2))',repelem(uG(1),length(unique(mzInd(uG==n,1:2))))'));
                end
            end
        end
    end
end
end

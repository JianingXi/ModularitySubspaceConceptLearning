function [TPR,FPR,AUC_ROC] = ROC_curve(TrueIDX,ScoreIDX)


    if size(TrueIDX,2) ~= 1
        TrueIDX = TrueIDX';
    end
    if size(ScoreIDX,2) ~= 1
        ScoreIDX = ScoreIDX';
    end
    
    [~,ind] = sort(ScoreIDX,'descend');
    RankOfTrue = TrueIDX(ind);
    FPR = cumsum(RankOfTrue == 0)/sum(RankOfTrue == 0);
    TPR = cumsum(RankOfTrue == 1)/sum(RankOfTrue == 1);
    AUC_ROC = sum(diff(FPR).*TPR(2:end));
end


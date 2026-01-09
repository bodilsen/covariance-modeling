function [RefreshData] = MakeRefresh2x2(Temp)

    [refPr,refTm,~] = realized_refresh_time_bivariate('seconds',Temp{1}(:,2),Temp{1}(:,1),Temp{2}(:,2),Temp{2}(:,1));
	RefreshData = [refTm refPr];
end

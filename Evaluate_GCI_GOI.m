function [IDR_GCI,MR_GCI,FAR_GCI,Bias_GCI,std_GCI,MSE_GCI,IDR_GOI,MR_GOI,FAR_GOI,Bias_GOI,std_GOI,MSE_GOI,VUDE,loc_fas_GCI,loc_missed_GCI] = Evaluate_GCI_GOI( sig,GCIs,GOIs,EstimatedGCIs,EstimatedGOIs,MinP,fs )
%EVALUATEGCI Summary of this function goes here
%   Detailed explanation goes here
[identsGCI, identsGOI, missesGCI, missesGOI, fasGCI, fasGOI, errorsGCI, ...
 errorsGOI, VUDE_samples,loc_fas_GCI,loc_missed_GCI] = Eval_GCI_GOI(GOIs,GCIs,EstimatedGOIs,...
                                         EstimatedGCIs,MinP,fs);

Total_NO_LarCyclesGCI = length(GCIs);
Total_NO_LarCyclesGOI = length(GOIs);

IDR_GCI = 100*(identsGCI./Total_NO_LarCyclesGCI);
IDR_GOI = 100*(identsGOI./Total_NO_LarCyclesGOI);
MR_GCI = 100*(missesGCI./Total_NO_LarCyclesGCI);
MR_GOI = 100*(missesGOI./Total_NO_LarCyclesGOI);
FAR_GCI = 100*(fasGCI./Total_NO_LarCyclesGCI);
FAR_GOI = 100*(fasGOI./Total_NO_LarCyclesGOI);
Bias_GCI = mean(errorsGCI);
Bias_GOI = mean(errorsGOI);
std_GCI  = std(errorsGCI);
std_GOI  = std(errorsGOI);
MSE_GCI = (1/length(errorsGCI))*sum(errorsGCI.^2);
MSE_GOI = (1/length(errorsGOI))*sum(errorsGOI.^2);
VUDE = 100*(VUDE_samples./length(sig));

end


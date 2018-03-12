function [mu,sig] = logn_dist_params(Mn,Mw)

mu=log(Mn)-1/2*log(Mw/Mn); sig = sqrt(log(Mw/Mn));

end



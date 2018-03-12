function [dp_dist, mw_dist] = generate_poly_dist(Mn,Mw,mon_weight)

mu=log(Mn)-1/2*log(Mw/Mn); sig = sqrt(log(Mw/Mn));

mw_dist = lognrnd(mu,sig,100000,1);
dp_dist = mw_dist./mon_weight;

end



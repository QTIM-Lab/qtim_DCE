function cost = kin_fit_test(params)    

FR_mins = .1;


t = 0:(50-1);
t = t.*FR_mins;

k1=exp(params(1));
    Ve=1/(1+exp(-params(2)));
    k2 = k1/Ve;
    
    res = exp(-k2.*t);
    est_conc = k1.*conv(AIF, res).*FR_mins;
    est_conc = est_conc(1:length(t));
    

    if size(est_conc,2)~=1
        est_conc=est_conc';
    end
 
    d = (obs_conc-est_conc).^2;
    cost = sum(d);
    
end
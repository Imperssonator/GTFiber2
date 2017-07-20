function fa_metrics = get_metrics_FA(fa_path,lcut)

fld_fa=FiberLengthsFA(fa_path);

load(fa_path)
% not_short = imageData.length_nm>lcut;
% imageData.xy = imageData.xy(not_short);
[op_fa,fibLenDens]=op2d_FA(imageData);
fa_metrics = [ ...
op_fa.Sfull;
op_fa.decayLen;
mean(fld_fa(fld_fa>lcut));
median(fld_fa(fld_fa>lcut));
fibLenDens];

end
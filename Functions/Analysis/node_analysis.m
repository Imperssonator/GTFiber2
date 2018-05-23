function [bpd, bp_im] = node_analysis(ims)

bpd = sum(sum(bwmorph(ims.skelTrim,'branchpoints'))) / ...
    prod(size(ims.gray).*ims.nmPix) * ...
    1000^2; % branch points per square micron

bp = bwmorph(ims.skelTrim,'branchpoints');
bp_im = ims.img;
bigbp = imdilate(bp,ones(5));
rr=bp_im(:,:,1);gg=bp_im(:,:,2);bb=bp_im(:,:,3);
rr(bigbp)=0;gg(bigbp)=255;bb(bigbp)=255;    % Cyan dots over nodes
bp_im=cat(3,rr,gg,bb);
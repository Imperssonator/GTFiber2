function skel2 = cleanSkel(skel,maxBranchLen)

B = bwmorph(skel, 'branchpoints');
E = bwmorph(skel, 'endpoints');
[m,n] = size(skel);
DB = bwdistgeodesic(skel,B);
DE = bwdistgeodesic(skel,E);
Nubs2 = DB+DE<maxBranchLen;
Nubs3 = Nubs2-B;
skel2 = skel-Nubs3;

NubsCheck = Nubs3 & DE>DB;

for ind = find(NubsCheck)'
    [i,j] = ind2sub([m,n],ind);
    nearBranchInd = findNearBranch(i,j,DB,m,n);
    if DE(i,j)>=DE(nearBranchInd)
        skel2(i,j)=1;
    end
end

skel2=bwareaopen(skel2,maxBranchLen);

end

function out = findNearBranch(i,j,DB,m,n)

if DB(i,j)==0 || isnan(DB(i,j)) || isinf(DB(i,j))
    out = sub2ind(size(DB),i,j);
    return
else
    hood = zeros(8,3);
    hood(:,1:2) = [i-1, j-1;...
            i-1, j;...
            i-1, j+1;...
            i, j-1;...
            i, j+1;...
            i+1, j-1;...
            i+1, j;...
            i+1, j+1];
    if i==1
        hood([1,2,3],:)=[];
    elseif i==m
        hood([6,7,8],:)=[];
    end
    if j==1
        hood([1,4,6],:)=[];
    elseif j==n
        hood([3,5,8],:)=[];
    end
    hood(:,3) = DB(sub2ind(size(DB),hood(:,1),hood(:,2)));
    [mind, minind] = min(hood(:,3),[],1);
    out = findNearBranch(hood(minind,1),hood(minind,2),DB,m,n);
end

end
function [L,C,allenergy] = kmeansplus(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 07/08/11
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = cumsum(sqrt(dot(D,D)));
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [tmp,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C).'));
    end
    iter=0;
    while any(L ~= L1) || iter==0
		iter=1;
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        %[tmp,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C).'),[],1);
		tmp=repmat(sum(X.^2,1),k,1)+repmat(sum(C.^2,1)',1,size(X,2))-2*C'*X;
		tmp2=cvx_cluster(tmp);
		for i=1:k
		L(find(tmp2(i,:)>=0.5,1))==i;
		end
	end
    
end
for i=1:k
    label{i}=find(L==i);
    temp=X(:,label{i})';
    energy(i)=sum(sum((temp-repmat(mean(temp),size(temp,1),1)).^2,2));
end
allenergy=sum(energy);
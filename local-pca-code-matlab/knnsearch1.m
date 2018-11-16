function idx=knnsearch1(point,X,radius_of_neighbor)
[N,D]=size(X);
X1=X-repmat(point,N,1);
dis=sum(X1.^2,2);
idx=find(dis<radius_of_neighbor^2);
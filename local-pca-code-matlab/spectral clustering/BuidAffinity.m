function A=BuidAffinity(X,option,sigma)
[N,D]=size(X);
normSquared=sum(X.^2,2);
disSquared=repmat(normSquared,1,N)+(repmat(normSquared,1,N))'-2*X*X';
switch option
case 1
A=exp(-disSquared/2/sigma^2);
case 2
scales=zeros(1,N);
for i=1:N
[a,b]=sort(disSquared(i,:),'ascend');
scales(i)=a(sigma);
end
scales=scales.^0.5;
A=exp(-disSquared./repmat(scales,N,1)./(repmat(scales,N,1))'/2);
case 3
scales=zeros(1,N);
for i=1:N
[a,b]=sort(disSquared(i,:),'ascend');
scales(i)=a(sigma);
end
scales=scales.^0.5;
temp=((2*scales'*scales)./(repmat(scales.^2,N,1)+(repmat(scales.^2,N,1))')).^0.5;
A=exp(-disSquared./(repmat(scales.^2,N,1)+(repmat(scales.^2,N,1))')).*temp;
end
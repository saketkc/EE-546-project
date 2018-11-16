function label=find_label(X,label1,mapping);
[N,D]=size(X);
label=zeros(1,N);
allnets=1:size(mapping,2);
%for i=1:N
dis=mapping(:,allnets);
[a,b]=min(dis');
label=label1(b);
%end

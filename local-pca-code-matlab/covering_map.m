function [ nets1 Dimss]=covering_map(X,neighborhood_size,neighborhood_type);
[U,S,V]=svd(X,0);
%X=U(:,1:4)*S(1:4,1:4);
[N,D]=size(X);
remain_labels=1:N;
switch neighborhood_type
case 5
center=X;
for i=1:N
nets{i}=knnsearch(center(i,:),X,neighborhood_size);
end
case 6
center=X;
for i=1:N
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
end
case 0
i=1;
center(i,:)=X(ceil(rand*N),:);
exclusion=knnsearch(center(i,:),X,neighborhood_size);
nets{i}=knnsearch(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
i=i+1;
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
exclusion=knnsearch(center(i,:),X,neighborhood_size);
nets{i}=knnsearch(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
end
case 1
i=1;
center(i,:)=X(ceil(rand*N),:);
exclusion=knnsearch1(center(i,:),X,neighborhood_size*2);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
i=i+1;
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
exclusion=knnsearch1(center(i,:),X,neighborhood_size*2);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
end
case 1.1
i=1;
center(i,:)=X(ceil(rand*N),:);
exclusion=knnsearch1(center(i,:),X,neighborhood_size);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
i=i+1;
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
exclusion=knnsearch1(center(i,:),X,neighborhood_size);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
end
case 1.15
i=1;
center(i,:)=X(ceil(rand*N),:);
exclusion=knnsearch(center(i,:),X,1);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
i=i+1;
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
exclusion=knnsearch(center(i,:),X,1);
nets{i}=knnsearch1(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
end
case 1.2
i=0;
center(1,:)=zeros(1,D);
while length(remain_labels)>0

center_temp=X(remain_labels(ceil(rand*length(remain_labels))),:);
[exclusion dim neighbors1]=knnsearch_adaptive(center_temp,X,0.1:0.005:0.25);
if dim==2
i=i+1;
center(i,:)=center_temp;
nets{i}=exclusion;remain_labels=setdiff(remain_labels,exclusion);
else
remain_labels=setdiff(remain_labels,exclusion);
end
end
case 2
i=1;
center(i,:)=X(ceil(rand*N),:);
exclusion=knnsearch(center(i,:),X,neighborhood_size*2);
nets{i}=knnsearch(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
i=i+1;
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
exclusion=knnsearch(center(i,:),X,neighborhood_size*2);
nets{i}=knnsearch(center(i,:),X,neighborhood_size);
remain_labels=setdiff(remain_labels,exclusion);
end
case 3
Y=zeros(size(X,1),neighborhood_size);
for i=1:size(X,1)
Y(i,:)=knnsearch(X(i,:),X,neighborhood_size);
end
i=1;
point=ceil(rand*N);
nets{i}=Y(point,:);
for j=nets{i}
remain_labels=setdiff(remain_labels,find((sum(Y'==j))==1));
end
while length(remain_labels)>0
i=i+1;
point=remain_labels(ceil(rand*length(remain_labels)));

nets{i}=Y(point,:);
for j=nets{i}
remain_labels=setdiff(remain_labels,find((sum(Y'==j))==1));
end
end
case 4
opts.MinNetPts=floor(0.01*N); 
opts.nPtsPerScale=floor(0.005*N);
opts.nScales=20;
opts.alpha0=0.2;
i=1;
center(i,:)=X(ceil(rand*N),:);
opts.seeds=center(i,:);
cd ../mapa
[noise estDims GoodScales radius(i)]=estimate_noise(X,opts.nPtsPerScale*opts.nScales+2*opts.MinNetPts,opts);
cd ../newproject
exclusion=knnsearch1(center(i,:),X,radius(i));
nets0{i}=exclusion;
remain_labels=setdiff(remain_labels,exclusion);
while length(remain_labels)>0
%if estDims<3
i=i+1;
%end
center(i,:)=X(remain_labels(ceil(rand*length(remain_labels))),:);
opts.seeds=center(i,:);

cd ../mapa
[noise estDims GoodScales radius(i)]=estimate_noise(X,opts.nPtsPerScale*opts.nScales+2*opts.MinNetPts,opts);
Dims(i)=estDims;
cd ../newproject
exclusion=knnsearch1(center(i,:),X,radius(i));
nets0{i}=exclusion;
remain_labels=setdiff(remain_labels,exclusion);
end
Dims;
nets=nets0;
% i=0;
% N0=size(center,1);
% remain_labels=1:N0;
% norm_center_squared=sum(center.^2,2);
% dis_center=real((repmat(norm_center_squared,1,N0)+(repmat(norm_center_squared,1,N0))'-2*center*center').^0.5);
% while length(remain_labels)>0
% i=i+1;
% [a,b]=max(radius(remain_labels));
% ind(i)=remain_labels(b);
% distance_table=dis_center(:,ind)-repmat(radius',1,i)-repmat(radius(ind),N0,1);
% if i>1
% alldistance=min(distance_table');
% else
% alldistance=distance_table';
% end
% remain_labels=find(alldistance>0);
% end
% ii=1;
% for i=1:length(ind)

% nets{ii}=nets0{ind(i)};
% ii=ii+1;
% Dimss(ii)=Dims(ind(i));

% end

for i=1:length(nets)
lengthnets(i)=length(nets{i});
end

j=1;
nets1{1}=1;
for i=1:length(nets)
	
    if lengthnets(i)>1
        nets1{j}=nets{i};j=j+1;
    end
end
end
if neighborhood_type~=4
nets1=nets;
Dimss=2*ones(1,length(nets1));
else
Dimss=Dims;
end



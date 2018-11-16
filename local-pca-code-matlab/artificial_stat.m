[N,D]=size(X);
X1=sum(X.^2,2);
dis=repmat(X1,1,N)+(repmat(X1,1,N))'-2*X*X';
R=(max(max(dis))).^0.5;

load twocurves0
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
twocurves0_stat(i)=allerror(4);
end

load twocurves1
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
twocurves1_stat(i)=allerror(4);
end

load twocurves2
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
twocurves2_stat(i)=allerror(4);
end

load twocurves3
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
twocurves3_stat(i)=allerror(4);
end

load twocurves4
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
twocurves4_stat(i)=allerror(4);
end

load threecurves
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.02,1.1);
threecurves_stat(i)=allerror(4);
end

load self_intersecting
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.1,1.1);
self_intersecting_stat(i)=allerror(4);
end

load twoeights
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,1,0.1,1.1);
twoeights_stat(i)=allerror(4);
end

load twospheres
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(X,M,2,0.2,1.1);
twospheres_stat(i)=allerror(4);
end

load toruscylinder
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(X,M,2,0.3,1.1);
toruscylinder_stat(i)=allerror(4);
end

load mobius
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,2,0.1,1.1); 
mobius_stat(i)=allerror(4);
end

load monkeysaddle
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Z,M,2,0.1,1.1); 
monkeysaddle_stat(i)=allerror(4);
end

load paraboloids
for i=1:100
[label coarse_label label1 models allerror allerror1 nets min_i min_j]=mani_clu(Y,M,2,0.07,1.1); 
paraboloids_stat(i)=allerror(4);
end


save artificial_stat twocurves1_stat twocurves2_stat twocurves3_stat twocurves4_stat threecurves_stat self_intersecting_stat twoeights_stat twospheres_stat toruscylinder_stat mobius_stat monkeysaddle_stat paraboloids_stat
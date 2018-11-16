function [label1 error] = SpectralClustering_main(X,M,option,sigma)
[N,D]=size(X);
 A=BuidAffinity(X,option,sigma);
 [label1] = SpectralClustering(A,max(M),1);
cd ..
for i=1:size(label1,1)
error(i)=testerror(label1(i,:),M)/length(M);
end
cd ./spectral' clustering'/
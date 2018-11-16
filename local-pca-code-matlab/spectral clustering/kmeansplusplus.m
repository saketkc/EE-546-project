function [allL,allC,allenergy] = kmeansplusplus(X,k,iter)
[allL,allC,allenergy] = kmeansplus(X,k);
for i=1:iter
    [L,C,energy] = kmeansplus(X,k);
    if energy<allenergy
        allL=L;
        allC=C;
        allenergy=energy;
    end
end    
    
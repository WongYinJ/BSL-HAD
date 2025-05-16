function [X] = CGra(E,F,Tau,X,B)
%CGRA 此处显示有关此函数的摘要
%   此处显示详细说明
R=B;
K=length(Tau)-1;
for k=1:K+1
    R=R-Tau(k)*E{k}*X*F{k};
end
P=R;
maxIter=1e3;
tol=1e-9;
for iter=1:maxIter
    aux=0;
    for k=1:K+1
        aux=aux+Tau(k)*trace(E{k}*P*F{k}*P');
    end
    alpha=R(:)'*R(:)/aux;
    X=X+alpha*P;
    aux=0;
    for k=1:K+1
        aux=aux+Tau(k)*E{k}*P*F{k};
    end
    R_new=R-alpha*aux;
    
    if norm(R_new)<tol
        break
    end
    beta=(R_new(:)'*R_new(:))/(R(:)'*R(:));
    P=R_new+beta*P;
    R=R_new;
end
[norm(R_new),iter];
end


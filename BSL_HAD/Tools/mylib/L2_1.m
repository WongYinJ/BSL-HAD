function [U] = L2_1(Y,lambda)
%L2_1 此处显示有关此函数的摘要
%   This function solves min_U: \lambda*||U||_{2,1}+1/2*||U-Y||_F^2
U=zeros(size(Y));
for j=1:size(U,2)
    U(:,j)=max(norm(Y(:,j),'fro')-lambda,0)/norm(Y(:,j),'fro')*Y(:,j);
end
end


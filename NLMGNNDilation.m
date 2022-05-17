function [output]=NLMGNNDilation(inputI,map)

[m,n]=size(inputI);
%R=kron(ones(m*n,1),inputI(:)');
% map(map>0)=1;
% RR=kron(ones(m*n,1),inputI(:)').*map;
MAPP1=map;
MAPP1(MAPP1>0)=1;
RR=kron(ones(m*n,1),inputI(:)').*MAPP1;
BB1=max(RR,[],2);
BB2=reshape(BB1,m,n);
output=full(BB2);%腐蚀输出，确实点由input补全

 
end

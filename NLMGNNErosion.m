function [output]=NLMGNNErosion(inputI,map)

inputI=255-inputI;
[output]=255 - NLMGNNDilation(inputI,map);


% [m,n]=size(inputI);

%R=kron(ones(m*n,1),inputI(:)');

% % map(map<=0)=1000;
% % RR=kron(ones(m*n,1),inputI(:)').*map;

% MAPP2=map;
% MAPP2(MAPP2<=0)=1000;
% 
% RR=kron(ones(m*n,1),inputI(:)').*MAPP2;
% 
% 
% BB1=min(RR,[],2);
% BB2=reshape(BB1,m,n);
% output=full(BB2);%腐蚀输出，确实点由input补全






% [m,n]=size(inputI);
% %R=kron(ones(m*n,1),inputI(:)');
% map(map>0)=1;
% R=inputI(:)';
% RR=kron(ones(m*n,1),R).*map;
% %r=sum(full(map(1,:)));
% RR=sort(RR,2);
% %sum(sum(RR))
% RRR=zeros(m*n,1);
% for i=1:m*n
%     waitbar(i/m/n);
%     temp=RR(i,:);
%     
%     
%     if sum(temp)>0
%         %find(temp,1)
%         RRR(i)=R(find(temp,1));
%         
%     end
%     
% end
% %sum(sum(R))
% 
% BB2=reshape(R,m,n);
% output=full(BB2);%腐蚀输出，确实点由input补全
 
end
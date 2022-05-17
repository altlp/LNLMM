function [B,W,WG]=fastKmeansfiltapprox2(A,S,h,Centre,Clusterindex,spatialkernel,convmethod,fast_flag,Aguide)
% Main filtering code lies here

if ~exist('Aguide','var')
     % guideimage is same as inputimage
     Aguide=A;
end
if ~exist('fast_flag','var')
     % guideimage is same as inputimage
     fast_flag=1;
end
guided=size(Aguide,3);
B=zeros(size(A));
Cluster=size(Centre,1);

%% Filtering using matlab command for convolutions

if strcmp(convmethod,'matlab')
    if strcmp(spatialkernel,'box')
        filt     = ones(2*S+1,2*S+1);       
    elseif strcmp(spatialkernel,'gaussian')       
        w  = round(6*S); if (mod(w,2) == 0); w  = w+1; end
        filt     = fspecial('gaussian', [w w], S);
    else
    end        
    for i=1:Cluster
        W=sum((Aguide-reshape(Centre(i,:),1,1,guided)).^2,3);      
        W=exp(-W./(2*h*h));        
        Wb=imfilter(W,filt);
        B=B+bsxfun(@times,bsxfun(@rdivide,imfilter(bsxfun(@times,W,A),filt),Wb),(Clusterindex==i));                 
    end
end    
%% Filtering using O(1) convolutions
m=size(A,1);
n=size(A,2);
CC=ones(m,n);
WG=zeros(m,n,Cluster);
if strcmp(convmethod,'O1') %nlm为01
    for i=1:Cluster
        W=sum((Aguide-reshape(Centre(i,:),1,1,guided)).^2,3);    %centre(i,;)代替搜索窗其他的块μ  Aguided为Apca
        W=exp(-W./(2*h*h));%φ函数
        WG(:,:,i)=W;
        if strcmp(spatialkernel,'box')%nlm中为box                                                                           %box_filter(bsxfun(@times,W,A)为
            Wb=box_filter(W,S,fast_flag);  %s=10    求ζ                                                                     %A为输入噪声图片   W*A为h(nlm未归一的分子项)
            B=B+bsxfun(@times,bsxfun(@rdivide,box_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(Clusterindex==i));
%             B=B+bsxfun(@times,bsxfun(@rdivide,box_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(1));
            % 对两个矩阵A和B之间的每一个元素进行指定的计算（函数fun指定） fast_flag=1 
        elseif strcmp(spatialkernel,'gaussian')
            Wb=gauss_filter(W,S,fast_flag);
            B=B+bsxfun(@times,bsxfun(@rdivide,gauss_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(Clusterindex==i));       
        else
        end 
    end
end 
end

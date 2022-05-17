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
if strcmp(convmethod,'O1') %nlmΪ01
    for i=1:Cluster
        W=sum((Aguide-reshape(Centre(i,:),1,1,guided)).^2,3);    %centre(i,;)���������������Ŀ��  AguidedΪApca
        W=exp(-W./(2*h*h));%�պ���
        WG(:,:,i)=W;
        if strcmp(spatialkernel,'box')%nlm��Ϊbox                                                                           %box_filter(bsxfun(@times,W,A)Ϊ
            Wb=box_filter(W,S,fast_flag);  %s=10    ���                                                                     %AΪ��������ͼƬ   W*AΪh(nlmδ��һ�ķ�����)
            B=B+bsxfun(@times,bsxfun(@rdivide,box_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(Clusterindex==i));
%             B=B+bsxfun(@times,bsxfun(@rdivide,box_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(1));
            % ����������A��B֮���ÿһ��Ԫ�ؽ���ָ���ļ��㣨����funָ���� fast_flag=1 
        elseif strcmp(spatialkernel,'gaussian')
            Wb=gauss_filter(W,S,fast_flag);
            B=B+bsxfun(@times,bsxfun(@rdivide,gauss_filter(bsxfun(@times,W,A),S,fast_flag),Wb),(Clusterindex==i));       
        else
        end 
    end
end 
end

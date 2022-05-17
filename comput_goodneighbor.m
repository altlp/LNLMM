WEI=WG2;
W2=WEI;%聚类后的存放的类别矩阵
[m,n,d2]=size(I);
t=10;
SS=100;%前20个邻居
IDXm=0; 
B=zeros(m*n*SS,5); 
tic
for i=1:m
     waitbar(i/m)
    for j=1:n
rmin = max(i-t,1);
rmax = min(i+t,m);
smin = max(j-t,1);
smax = min(j+t,n);
CC=(j-1)*m+i;
dd=Clusterindex(CC);
% ff=W2(:,:,dd);%这种赋值方式导致速度很慢
% ff=W2(:,:,1);
A=W2(rmin:1:rmax,smin:1:smax,dd);
% A=ff(rmin:1:rmax,smin:1:smax);
A=A./sum(sum(A));
A2=A';
[x2,x1]=meshgrid(rmin:1:rmax,smin:1:smax);
D=[x2(:),x1(:),A2(:)];
D=sortrows(D,3,'descend');
SSD=min(size(D,1),SS);
D=D(1:SSD,:);
fd=sum(D(:,3));
 D(:,3)=D(:,3)./fd;
B(IDXm+1:IDXm+SSD,:)=[i*ones(SSD,1),j*ones(SSD,1), D];  
IDXm=IDXm+SSD;      
    end  
end
toc
B=B([1:IDXm],:);

[m,n,d2]=size(I);
tic
U=(B(:,2)-1)*m+B(:,1);
V=(B(:,4)-1)*m+B(:,3); %
nlmedge=[U,V];
SW=B(:,5);

edgesP=[U,V];
N=max(max(edgesP));
% MAP0=adjacency(edgesP,SW,N);
MAP0=sparse(U,V,SW,N,N);%带权重的稀疏矩阵
MAP1=sqrt(MAP0 .* MAP0');%好邻居邻接矩阵图，得到的好邻居
MAP=MAP0.*(MAP1>0);%好邻居权值矩阵图
I3=I2(:,:,1);
output0=MAP*I3(:);
MAP2=sum(MAP,2);
output1=output0./MAP2;
% point0=find(MAP2==0);
% output1(point0)=Iact(point0)./255;
output=reshape(output1,[m,n]);
figure,imshow(output,[])
toc


function  [phs,phs_wrapped,ah] = PUNF(II)

image_size=size(II);
dimension=numel(image_size);
if dimension==3
    II = rgb2gray(II);
end
II=double(II);    


ff=fftshift(fft2(II));
FF=abs(fftshift(fft2(II)));           %计算全息图的频谱

% figure,imshow(FF,[0,max(max(FF))/100],'border','tight');%显示计算全息图的频谱
[r,c] = size(ff);

mid_r = floor(r / 2)+1;
mid_c = floor(c / 2)+1;

h = zeros(r,c); 


ra = 599;

%d=55 for RBC, d=87 for PS microsphere
d = 55;
%d = 87;



ra_r = mid_r - 3;
ra_c = mid_c + ra; 


max_value = 0;

max_r = 0;
max_c = 0;


%找出最大值和相应坐标
for i = ra_r - d  : ra_r + d  
    for j = ra_c - d : ra_c  + d
        if(FF(i,j) >= max_value)
            max_value = FF(i,j);
%             max_all = [max_all,max];
%             max_r = [max_r,i];
%             max_c = [max_c,j];
            max_r = i;
            max_c = j;
        end
    end
end



for i = max_r - d  : max_r + d  
    for j = max_c - d : max_c  + d
        if((i - max_r)^2 + (j - max_c)^2 <= d^2)
            h(i,j) =1;
        end
    end
end


h = ff.* h;
% ha = abs(h);
% figure,imshow(ha,[0,max(max(ha))/100],'border','tight');%显示计算全息图的频谱


h_tmp = zeros(r,c);
h_tmp(mid_r - d : mid_r + d,mid_c - d : mid_c + d) = h(max_r - d : max_r + d ,max_c - d : max_c + d);

h = h_tmp;

h_if = ifft2(ifftshift(h));
ah = sqrt(h_if.*conj(h_if));
ph=angle(h_if);   
phs_wrapped = ph;

a=ph;                                             %将包裹相位ph赋值给a
%下面开始进行最小二乘解包裹运算
[M,N]=size(a);                                    %计算二维包裹相位的大小（行、列数）
dx=zeros(M,N);dy=zeros(M,N);                      %预设包裹相位沿x方向和y方向的梯度
m=1:M-1; 
dx(m,:)=a(m+1,:)-a(m,:);                          %计算包裹相位沿x方向的梯度
dx=dx-pi*round(dx/pi);                            %去除梯度中的跳跃
n=1:N-1;
dy(:,n)=a(:,n+1)-a(:,n);                          %计算包裹相位沿y方向的梯度
dy=dy-pi*round(dy/pi);                            %去除梯度中的跳跃
p=zeros(M,N);p1=zeros(M,N);p2=zeros(M,N); %为计算ρnm作准备
m=2:M;
p1(m,:)=dx(m,:)-dx(m-1,:);                        %计算Δgxnm-Δgx(n-1)m
n=2:N;
p2(:,n)=dy(:,n)-dy(:,n-1);                        %计算Δgynm–Δgyn(m-1)
p=p1+p2;                                          %计算ρnm
p(1,1)=dx(1,1)+dy(1,1);                           %计算ρnm
n=2:N;
p(1,n)=dx(1,n)+dy(1,n)-dy(1,n-1);                 %赋值Neumann边界条件
m=2:M;
p(m,1)=dx(m,1)-dx(m-1,1)+dy(m,1);
pp=dct2(p)+eps;                                   %计算ρnm的DCT
fi=zeros(M,N);
for m=1:M                                         %计算Φnm在DCT域的精确解
   for n=1:N  
      fi(m,n)=pp(m,n)/(2*cos(pi*(m-1)/M)+2*cos(pi*(n-1)/N)-4+eps);
   end
end
fi(1,1)=pp(1,1);                                  %赋值DCT域的Φ11
phs=idct2(fi);                                    %用iDCT计算解包裹相位在空域中的值



end


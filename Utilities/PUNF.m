function  [phs,phs_wrapped,ah] = PUNF(II)

image_size=size(II);
dimension=numel(image_size);
if dimension==3
    II = rgb2gray(II);
end
II=double(II);    


ff=fftshift(fft2(II));
FF=abs(fftshift(fft2(II)));           %����ȫϢͼ��Ƶ��

% figure,imshow(FF,[0,max(max(FF))/100],'border','tight');%��ʾ����ȫϢͼ��Ƶ��
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


%�ҳ����ֵ����Ӧ����
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
% figure,imshow(ha,[0,max(max(ha))/100],'border','tight');%��ʾ����ȫϢͼ��Ƶ��


h_tmp = zeros(r,c);
h_tmp(mid_r - d : mid_r + d,mid_c - d : mid_c + d) = h(max_r - d : max_r + d ,max_c - d : max_c + d);

h = h_tmp;

h_if = ifft2(ifftshift(h));
ah = sqrt(h_if.*conj(h_if));
ph=angle(h_if);   
phs_wrapped = ph;

a=ph;                                             %��������λph��ֵ��a
%���濪ʼ������С���˽��������
[M,N]=size(a);                                    %�����ά������λ�Ĵ�С���С�������
dx=zeros(M,N);dy=zeros(M,N);                      %Ԥ�������λ��x�����y������ݶ�
m=1:M-1; 
dx(m,:)=a(m+1,:)-a(m,:);                          %���������λ��x������ݶ�
dx=dx-pi*round(dx/pi);                            %ȥ���ݶ��е���Ծ
n=1:N-1;
dy(:,n)=a(:,n+1)-a(:,n);                          %���������λ��y������ݶ�
dy=dy-pi*round(dy/pi);                            %ȥ���ݶ��е���Ծ
p=zeros(M,N);p1=zeros(M,N);p2=zeros(M,N); %Ϊ�����nm��׼��
m=2:M;
p1(m,:)=dx(m,:)-dx(m-1,:);                        %���㦤gxnm-��gx(n-1)m
n=2:N;
p2(:,n)=dy(:,n)-dy(:,n-1);                        %���㦤gynm�C��gyn(m-1)
p=p1+p2;                                          %�����nm
p(1,1)=dx(1,1)+dy(1,1);                           %�����nm
n=2:N;
p(1,n)=dx(1,n)+dy(1,n)-dy(1,n-1);                 %��ֵNeumann�߽�����
m=2:M;
p(m,1)=dx(m,1)-dx(m-1,1)+dy(m,1);
pp=dct2(p)+eps;                                   %�����nm��DCT
fi=zeros(M,N);
for m=1:M                                         %���㦵nm��DCT��ľ�ȷ��
   for n=1:N  
      fi(m,n)=pp(m,n)/(2*cos(pi*(m-1)/M)+2*cos(pi*(n-1)/N)-4+eps);
   end
end
fi(1,1)=pp(1,1);                                  %��ֵDCT��Ħ�11
phs=idct2(fi);                                    %��iDCT����������λ�ڿ����е�ֵ



end


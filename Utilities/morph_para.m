function data = morph_para(d_sam,proc_mask)
    stats_new = regionprops(proc_mask, 'BoundingBox' ,'Area','Centroid' ,'PixelList' ,'Perimeter','Eccentricity','Circularity');
    [~,num] = bwlabel(proc_mask,8);

    PA_ave = zeros(num,1);%用于记录每个细胞的投影面积
    height_ave = zeros(num,1); %用于记录每个细胞的平均厚度
    volume_ave = zeros(num,1);%用于记录每个细胞的体积
    surface_ave = zeros(num,1);%用于记录每个细胞的表面积
    sphericity_ave = zeros(num,1);%用于记录每个细胞的球形度
    MCD_ave = zeros(num,1);%用于记录每个细胞的最小圆柱直径
    Drymass_ave = zeros(num,1);%用于记录每个细胞的干质量
    kurtosis_ave = zeros(num,1);%用于记录每个细胞的高度峰度
    skewness_ave = zeros(num,1);%用于记录每个细胞的高度偏度
    perimeter_ave = zeros(num,1);%用于记录每个细胞的周长
    diameter_ave = zeros(num,1);%用于记录每个细胞的圆直径
    eccentricity_ave = zeros(num,1);%用于记录每个细胞的偏心率
    circularity_ave = zeros(num,1);%用于记录每个细胞的圆度
    data = zeros(num,13); %13个参数

    pixel = 2.5/40;% 像素宽度为2.5um，但是放大了40倍，所以这里除以40倍方便后续求细胞实际大小

    for i=1:num
        area = stats_new(i).Area; %连通区域的面积
        pointList = stats_new(i).PixelList; %每个连通区域的像素位置
        rIndex=pointList(:,2);cIndex=pointList(:,1);
        labelcell = zeros(length(rIndex),1);   %用于记录某一连通区域对应的所有像素的值
        for s = 1:length(rIndex)
            labelcell(s,1) = d_sam(rIndex(s),cIndex(s));
        end
        labelcell = labelcell - min(min(labelcell));%把该区域对应的高度值分布设置为从0开始
        height_ave(i) = mean(labelcell);
        PA_ave(i) = area*pixel*pixel;
        volume_ave(i) = pixel*pixel*sum(labelcell(:));
        surface_ave(i) = PA_ave(i)+mySurfaceArea(rIndex,cIndex,d_sam,pixel);
        sphericity_ave(i) = 4*pi*(3*volume_ave(i)/(4*pi))^(2/3)/surface_ave(i);
        p = pi/12;
        q = 0;
        c = -surface_ave(i)/4;
        d = volume_ave(i);
        MCD_ave(i) = Solve3Polynomial(p,q,c,d);
        alpha = 0.2*1e6;
        delta_n = 0.06;
        Drymass_ave(i) = area*pixel*pixel/alpha*delta_n*mean(labelcell);
        kurtosis_ave(i) = kurtosis(labelcell);
        skewness_ave(i) = skewness(labelcell);
        perimeter_ave(i) = pixel*stats_new(i).Perimeter;
        diameter_ave(i) = 2*sqrt(PA_ave(i)/pi);
        eccentricity_ave(i) = stats_new(i).Eccentricity;
        circularity_ave(i) = stats_new(i).Circularity;
    end
    data = [height_ave PA_ave volume_ave surface_ave sphericity_ave MCD_ave Drymass_ave kurtosis_ave skewness_ave perimeter_ave diameter_ave eccentricity_ave circularity_ave];

end


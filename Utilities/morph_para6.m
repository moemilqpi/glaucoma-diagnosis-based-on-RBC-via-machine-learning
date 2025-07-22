function data = morph_para6(d_sam,proc_mask, visual)
    stats_new = regionprops(proc_mask, 'BoundingBox' ,'Area','Centroid' ,'PixelList' );
    [~,num] = bwlabel(proc_mask,8);

    PA_ave = zeros(num,1);%用于记录每个细胞的投影面积
    height_ave = zeros(num,1); %用于记录每个细胞的平均厚度
    volume_ave = zeros(num,1);%用于记录每个细胞的体积
    surface_ave = zeros(num,1);%用于记录每个细胞的表面积
    sphericity_ave = zeros(num,1);%用于记录每个细胞的球度
    MCD_ave = zeros(num,1);%用于记录每个细胞的最小圆柱直径
    Drymass_ave = zeros(num,1);%用于记录每个细胞的干质量
    kurtosis_ave = zeros(num,1);%用于记录每个细胞的峰度
    skewness_ave = zeros(num,1);%用于记录每个细胞的偏度
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
        % p = -surface_ave(i)*3/pi;
        % q = volume_ave(i)*12/pi;
        % MCD_ave(i) = (-q/2+sqrt((q/2)^2+(p/3)^3))^(1/3)+(-q/2-sqrt((q/2)^2+(p/3)^3))^(1/3);
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
    end

    MeanHeightImage = mean(height_ave);
    StdHeightImage = std(height_ave);
    MeanPAImage = mean(PA_ave);
    StdPAIImage = std(PA_ave);
    MeanVolumeImage = mean(volume_ave);
    StdVolumeImage = std(volume_ave);
    MeanSurfaceImage = mean(surface_ave);
    StdSurfaceImage = std(surface_ave);
    MeanSphericityImage = mean(sphericity_ave);
    StdSphericityImage = std(sphericity_ave);
    MeanMCDImage = mean(MCD_ave);
    StdMCDImage = std(MCD_ave);
    MeanDrymass = mean(Drymass_ave);
    StdDrymass = std(Drymass_ave);
    MeanKurtosis = mean(kurtosis_ave);
    StdKurtosis = std(kurtosis_ave);
    MeanSkewness = mean(skewness_ave);
    StdSkewness = std(skewness_ave);
    display(num);
    display(MeanHeightImage);
    display(MeanPAImage);
    display(MeanVolumeImage);
    display(MeanSurfaceImage);
    display(MeanSphericityImage);
    display(MeanMCDImage);
    display(MeanDrymass);
    display(MeanKurtosis);
    display(MeanSkewness);
    
    %将数据导入到excel中
    data = {visual, num, MeanHeightImage, StdHeightImage, MeanPAImage, StdPAIImage, MeanVolumeImage, StdVolumeImage,...
    MeanSurfaceImage, StdSurfaceImage, MeanSphericityImage, StdSphericityImage, MeanMCDImage, StdMCDImage,...
    MeanDrymass, StdDrymass, MeanKurtosis, StdKurtosis, MeanSkewness, StdSkewness};

end


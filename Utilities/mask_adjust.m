function [mask_new] = mask_adjust(phase,mask,threshold,min_projected_area,max_projected_area,phase_range)
    [~,num] = bwlabel(mask,8);
    [B,~] = bwboundaries(mask,'noholes');
    stats = regionprops(mask, 'BoundingBox' ,'Area','Centroid' ,'PixelList' ); %统计白色的连通区域
    
    for i = 1:num
        area = stats(i).Area; % The area of the connected region
        pointList = stats(i).PixelList; %每个连通区域的像素位置
        rIndex=pointList(:,2);cIndex=pointList(:,1);
        pointLength = length(rIndex);
        
        %圆形度区域检索
          % obtain (X,Y) boundary coordinates corresponding to label 'k'
          boundary = B{i};
    
          % compute a simple estimate of the object's perimeter
          delta_sq = diff(boundary).^2;    
          perimeter = sum(sqrt(sum(delta_sq,2)));
    
          % compute the roundness metric
          metric = 4*pi*area/perimeter^2;
          % 相位值差值计算
            phase_values = zeros(pointLength, 1);
            for k = 1:pointLength
                phase_values(k) = phase(rIndex(k), cIndex(k));
            end
            phase_dif = max(phase_values) - min(phase_values);
          
          % 圆心和半径估计
          centroid = stats(i).Centroid; % [x, y]
          distances = sqrt((cIndex - centroid(1)).^2 + (rIndex - centroid(2)).^2);
          radius_est = max(distances);
          outer_indices = find(distances == radius_est / 2);
 
        if isempty(outer_indices)
            avg_outer_phase = mean(phase_values); % fallback
        else
            avg_outer_phase = mean(phase_values(outer_indices));
        end

        center_phase = phase(round(centroid(2)), round(centroid(1)));
        
        if (area < min_projected_area  || area > max_projected_area || metric < threshold || phase_dif < phase_range || avg_outer_phase <= center_phase)
            for x = 1:pointLength
                    mask(rIndex(x),cIndex(x)) = 0;
            end
        end
    end

    mask_new = mask;
    
end


clc
clear all
close all

%% after the segmentation is completed using the SAM2 model, further processing of the mask is carried out, including setting the cell area threshold, setting the roundness, etc. 
% Save the processed mask and superimpose it on the surface of the phase diagram to display the effect
glaucoma_num = 26;
healthy_num = 26;
volunteer_num = glaucoma_num + healthy_num;
img_width = 2160;
img_height = 2600;
for t = 1:volunteer_num   
    if (t<healthy_num+1)  
        volunteer = 'Healthy\';
        mask_folder = ['Data\', volunteer, 'No.', num2str(t)];
    else
        volunteer = 'Glaucoma\';
        mask_folder = ['Data\', volunteer, 'No.', num2str(t-healthy_num)];
    end
    proc_mask_folder = [mask_folder,'\Proc_masks'];
    mkdir(proc_mask_folder);
    proc_mask_drawing_folder = [mask_folder, '\Proc_mask_drawing'];
    mkdir(proc_mask_drawing_folder);

    files = dir(strcat(mask_folder,'\masks','\*.png'));
    num = numel(files);
    for i=1:num
        %  Find the problematic phase. If meets the encounter one, skip it directly and do not handle the mask
        phase_folder = [mask_folder, '\phase maps', '\phase', num2str(i)];
        phase = struct2array(load(phase_folder));
        mean_left = mean(phase(:,1));
        mean_right = mean(phase(:,img_height));
        mean_up = mean(phase(1,:));
        mean_down = mean(phase(img_width,:));
        if (abs(mean_left-mean_right)>4) || (abs(mean_up-mean_down)>4)
            continue;
        end
        
        single_mask_folder = [mask_folder, '\masks', '\mask', num2str(i), '.png'];
        mask = imread(single_mask_folder);
        mask = logical(mask);
        circle_threshold = 0.8;% The minimum value of the circularity threshold
        min_area_threshold = 6000; % The minimum number of projected pixels
        max_area_threshold = 28000;% The maximum number of projected pixels
        phase_range = 0.7; % The phase range from the middle to the edge of the cell
        proc_mask = mask_adjust(phase,mask,circle_threshold,min_area_threshold, max_area_threshold, phase_range);
        
        
        proc_mask_name = [proc_mask_folder, '\proc_mask',num2str(i),'.png'];
        imwrite(uint8(proc_mask), proc_mask_name);
        
        %Draw the contour line
        proc_mask = logical(proc_mask);
        image_name = [mask_folder, '\images', '\image', num2str(i), '.png'];
        original_img = imread(image_name);
        if size(original_img, 3) == 1
            original_img = repmat(original_img, [1 1 3]);
        end
        outline = bwperim(proc_mask);

        se = strel('disk', 1);  
        thick_outline = imdilate(outline, se); 
        outline_color = [255, 0, 0];
        overlay_img = original_img;
        [row, col] = find(thick_outline);          
        for d = 1:length(row)
            overlay_img(row(d), col(d), :) = reshape(outline_color, [1 1 3]);
        end
        proc_mask_drawing_name = [proc_mask_drawing_folder, '\proc_mask_drawing', num2str(i), '.png'];
        imwrite(overlay_img, proc_mask_drawing_name);
        
    end
end
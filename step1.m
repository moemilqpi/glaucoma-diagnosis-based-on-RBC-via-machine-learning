clc
clear all
close all

%% The first step is to reconstruct all phase maps. The original phase maps are saved as mat files, and the normalized phase maps are saved in png format for the segmentation of the SAM2 model
glaucoma_num = 26;
healthy_num = 26;
volunteer_num = glaucoma_num + healthy_num;
for t = 1:volunteer_num
    if (t<healthy_num+1)
        volunteer = 'Healthy\';
        name1 = ['Data\',volunteer,'No.', num2str(t)];
        newfolder_img = [name1, '\images'];
        mkdir(newfolder_img);
        newfolder_phase = [name1, '\phase maps'];
        mkdir(newfolder_phase);
        newfolder_mask = [name1, '\masks'];
        mkdir(newfolder_mask);
    else
        volunteer = 'Glaucoma\';
        name1 = ['Data\',volunteer,'No.', num2str(t-healthy_num)];
        newfolder_img = [name1, '\images']; 
        mkdir(newfolder_img);
        newfolder_phase = [name1, '\phase maps'];
        mkdir(newfolder_phase);
        newfolder_mask = [name1, '\masks'];
        mkdir(newfolder_mask);
    end

    name_inter = '\inter';
    name_back = '\back';
    name_bmp = '.bmp';
    
    files = dir(strcat(name1,'\*.bmp'));
    num = numel(files);


    for i=1:num/2
        path_inter = strcat(name1,name_inter,num2str(i),name_bmp);
        path_back = strcat(name1,name_back,num2str(i),name_bmp);
        II_inter = imread(path_inter);
        II_back = imread(path_back);
        phs1 = PUNF(II_inter);
        phs2  = PUNF(II_back);   
        phs_last = phs1 - phs2;
        normalized = mat2gray(phs_last);  % Normalize phase_last to [0,1]

        name_phase = [newfolder_phase,'\', 'phase', num2str(i), '.mat'];
        save(name_phase,'phs_last');
        name_forseg = [newfolder_img,'\', 'image', num2str(i), '.png'];
        imwrite(normalized, name_forseg);
    end
end
%% Choose the slice and the area
figure; 

for i = 1:50

    slice  = i;
    area_x = 5:25; 
    area_y = 43:62;

%     aux = abs(u{1,1}{slice}(area_x,area_y));
    aux = abs(u_clean{slice}(area_x,area_y));
%     aux_gt = fMRI_groundTruth{slice + 10}(area_x,area_y);

    % Normalize 
    aux = aux/max(aux(:));
    
    surf(aux); 
    view(130,48); caxis([0,1]);
    
   % set(gca, 'visible', 'off')
    set(gcf, 'color', 'w');
    export_fig(['images/area_gif/gt_',num2str(slice)]);
    
end

%%
for j = 1:5

    for i = 1:10
    slice  = i;
    area_x = 5:25; 
    area_y = 43:62;

    aux = abs(u{1,j}{slice}(area_x,area_y));
%     aux_gt = fMRI_groundTruth{slice + 10}(area_x,area_y);

    % Normalize 
    aux = aux/max(aux(:));
    
    surf(aux); 
    view(130,48); caxis([0,1]);
    
   % set(gca, 'visible', 'off')
    set(gcf, 'color', 'w');
    export_fig(['images/area_gif/recon_',num2str(slice)]);
    end
    
end




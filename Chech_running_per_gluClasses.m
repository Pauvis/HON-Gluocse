% Find location of clustered cells
%%%%%% Analysis workflow:
% cell_f_extraction_2P_doublecheckplanes_24012022a.m
% ExportPreprocesseddata_2PGlusensstim08022022.m
% Combine_2Pdata0802.m
% Analysis_combined_data_2P.m
% Locate2Pcells.m
% Chech_running_per_gluClasses.m
clear all
close all

clc
% params.colour_list = {[0 0.4470 0.7410],[ 0.8500 0.3250 0.0980],[ 0.9290 0.6940 0.1250],[ 0.4940 0.1840 0.5560],[ 0.4660 0.6740 0.1880],[ 0.3010 0.7450 0.9330],[ 0.6350 0.0780 0.1840],[0,0,0],[0.5,0,0],[0,0.5,0],[0,0,0.5]};
colorlist = [0,0,0;100,143,255;120,94,240;220,38,127;254,97,0;255,176,0]/255;
colorlist = [0,13,255;40,120,255;255,30,30;255,100,0;0,0,0;70,62,70;]/255;
colorlist = [10,33,210;86,180,255;216,37,37;255,120,100;0,0,0;70,62,70;]/255; %
colorlist = [255,30,30;10,33,210;10,10,10;86,180,255;10,10,10;0,0,0;70,62,70;]/255;
colorlist = [216,37,37;255,120,100;10,33,210;86,180,255;10,10,10;0,0,0;70,62,70;]/255;% colorlist = [10,33,210;;10,10,10;216,37,37;255,120,100;0,0,0;70,62,70;]/255;

dir_t = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp';
fid = 'All_data.mat';
animal_list = {'077','082'}; % rotate by -45*
animal_list2 = {'126b','152f','G35','G36','m143a'}; %rotate by 180*

load(fullfile(dir_t,fid)); % main data file
%% add metadata to the experiments 
% For each animal
meta_data_folder = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp\Recording meta data\';
all_available_metafiles = dir(fullfile(meta_data_folder,'*.ini'));
all_available_metafiles = {all_available_metafiles.name};
for animal_i = 1:height(All_data2)
    % for each condition
    conditions = All_data2.Properties.VariableNames(2:end);
    for cond_i = 1:length(conditions)
        % find experiment names and itterate
        exp_names = fieldnames(All_data2.(conditions{cond_i}){animal_i});
        for exp_i = 1:length(exp_names)
            All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i});
            % properties to look for
            ID_t = All_data2.ID{animal_i};
            ID_t = strrep(ID_t,'m','');
            treatment = conditions{cond_i};
            exp_date = exp_names{exp_i}(5:12);
            
            % find corresponding meta file if available
            idx = find(contains(all_available_metafiles,ID_t,'IgnoreCase',true) & contains(all_available_metafiles,treatment,'IgnoreCase',true) & contains(all_available_metafiles,exp_date));
            if ~isempty(idx) && length(idx) == 1
                fid=fopen(fullfile(meta_data_folder,all_available_metafiles{idx}));
                t=textscan(fid,'%s','delimiter',' ');
                fid=fclose(fid);
                t = t{:};

                All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i}).zoom = str2num(t{find(strcmpi(t,'zoom'))+2});
                All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i}).noofplanes = str2num(t{find(strcmpi(t,'no.of.planes'))+2});
                All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i}).zdistance = str2num(t{find(strcmpi(t,'total.z.distance'))+2});
                All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i}).zspacing = str2num(t{find(strcmpi(t,'z.spacing'))+2});
                All_data2.(conditions{cond_i}){animal_i}.(exp_names{exp_i}).gridsize = str2num(strrep(strrep(t{find(strcmpi(t,'grid.size'))+2},'u',''),'"',''));
                
                
            elseif length(idx) > 1
                disp(strjoin({'Exp match 2 files for',ID_t,treatment,exp_date}));
                adads
            else
                disp(strjoin({'Exp details not found for',ID_t,treatment,exp_date}));
            end
        end
    end
    
end
%% adds cell classification based on glu responses
classes_originally = {'G', 'iG', 'dG', 'idG', 'None'}; % [0,1,2,3,4] <- this is how they were classified by Alex
classes = {'G', 'dG', 'iG', 'idG', 'None'}; % [0,1,2,3,4] <- This is the order in which they will be plotted
[~,~,IB] = intersect(classes_originally,classes,'stable');

classes_fid = 'Classifications.xlsx';
dir_tt = dir_t;
% new cell classification
dir_t2 = 'P:\Alexander\GlucoseAnalysis2_0\Figs\April2022';
dir_t2 = 'P:\Alexander\GlucoseAnalysis2_0\Figs\May2022';
classes_fid2 = 'ClassificationsNew.xlsx';

sheets = sheetnames(fullfile(dir_tt,classes_fid));

groups = []; groupsnos =[];
coordinate_map = []; coordinate_map_animal  = {}; groups_animal = {}; all_cells_map = {};
varnames = fieldnames(All_data2);
varnames = varnames(2:3);
for sheet_i = 1:length(sheets)
    classes_mat = xlsread(fullfile(dir_t,classes_fid),sheets{sheet_i});
    
    classes_mat2 = xlsread(fullfile(dir_t2,classes_fid2),sheets{sheet_i});
    % find experiments with info
    which_animal = find(strcmpi(All_data2.ID,sheets{sheet_i}) == 1);
    
    [~,sort_cells] = sort(classes_mat(:,1));
    classes_mat = classes_mat(sort_cells,:);
    classes_mat = classes_mat+1; % Fix python indexing
    
    [~,sort_cells] = sort(classes_mat2(:,1));
    classes_mat2 = classes_mat2(sort_cells,:);
    classes_mat2 = classes_mat2+1; % Fix python indexing
    
    % Remap based on new class order
    idxes_to_map = nan(size(classes_mat,1),length(IB));
    for i = 1:length(IB)
        idxes_to_map(:,i) = (classes_mat(:,2) == i);
    end
    for i = 1:length(IB)
        classes_mat(idxes_to_map(:,i)==1,2) = IB(i);
    end
    
    idxes_to_map = nan(size(classes_mat2,1),length(IB));
    for i = 1:length(IB)
        idxes_to_map(:,i) = (classes_mat2(:,2) == i);
    end
    for i = 1:length(IB)
        classes_mat2(idxes_to_map(:,i)==1,2) = IB(i);
    end
    
    
    exp_names = fieldnames(All_data2.Glu{which_animal}); %use the first exp only? mean?
    exp_names_sal = fieldnames(All_data2.Sal{which_animal}); %use the first exp only? mean?
    clear classes_tbl
    classes_tbl.cell_IDs = classes_mat(:,1);
    classes_tbl.cell_class = classes_mat(:,2);
    
    for class_i = 1:length(unique(classes_tbl.cell_class))
        idx = find(classes_tbl.cell_class == class_i);
        
        classes_tbl.cell_class_name(idx,1) = classes(class_i);
%         classes_tbl.classes =
    end
%     %%
    classes_tbl = struct2table(classes_tbl);
    % for each cell
    for cell_i = 1:size(classes_mat,1)
        
        classes_tbl.X{cell_i} = [];
        classes_tbl.Y{cell_i} = [];
        classes_tbl.Zoom{cell_i} = [];
        classes_tbl.Z{cell_i} = [];
        % find IDs.
%         classes_tbl.cell_identities{cell_i} = All_data2.Glu{which_animal}.(exp_names{1}).cell_identities{cell_i};
        % find section number
        classes_tbl.sections(cell_i) = All_data2.Glu{which_animal}.(exp_names{1}).plane_list(cell_i);
        % find location map
        for exp_i = 1:length(exp_names)
            classes_tbl.maps{cell_i,exp_i} = All_data2.Glu{which_animal}.(exp_names{exp_i}).ROI_mat{cell_i};
            classes_tbl.X{cell_i} = [classes_tbl.X{cell_i}; sum(All_data2.Glu{which_animal}.(exp_names{1}).ROI_mat{cell_i},1)];
            classes_tbl.Y{cell_i} = [classes_tbl.Y{cell_i}; sum(All_data2.Glu{which_animal}.(exp_names{1}).ROI_mat{cell_i},2)'];
            classes_tbl.Zoom{cell_i}    = [classes_tbl.Zoom{cell_i};    All_data2.Glu{which_animal}.(exp_names{1}).zoom];
            classes_tbl.Z{cell_i}       = [classes_tbl.Z{cell_i};       All_data2.Glu{which_animal}.(exp_names{1}).zdistance];
        end
        for exp_i = 1:length(exp_names_sal)
            classes_tbl.maps{cell_i,length(exp_names)+exp_i} = All_data2.Sal{which_animal}.(exp_names_sal{exp_i}).ROI_mat{cell_i};
            classes_tbl.X{cell_i} = [classes_tbl.X{cell_i}; sum(All_data2.Sal{which_animal}.(exp_names_sal{1}).ROI_mat{cell_i},1)];
            classes_tbl.Y{cell_i} = [classes_tbl.Y{cell_i}; sum(All_data2.Sal{which_animal}.(exp_names_sal{1}).ROI_mat{cell_i},2)'];
            classes_tbl.Zoom{cell_i}    = [classes_tbl.Zoom{cell_i};    All_data2.Sal{which_animal}.(exp_names_sal{1}).zoom];
            classes_tbl.Z{cell_i}       = [classes_tbl.Z{cell_i};       All_data2.Sal{which_animal}.(exp_names_sal{1}).zdistance];
        end
        classes_tbl.Zoom{cell_i} = median(classes_tbl.Zoom{cell_i} );
        classes_tbl.Z{cell_i} = median(classes_tbl.Z{cell_i} );

    end

%     %% combine ROIs into cells
    % combine cell lists for a given animal from all of the experiments
    cell_ID_mat = {};
    cell_IDs_comb = {};
    cell_f_std = [];
    cell_f_max = [];
    var_exp_used = {};
    for var_i = 1:length(varnames)
        
        exp_list = fieldnames(All_data2.(varnames{var_i}){sheet_i});
        plane_list = All_data2.(varnames{var_i}){sheet_i}.(exp_list{1}).plane_list;
        count = 1;
        for exp_i = 1:length(exp_list)

            t = All_data2.(varnames{var_i}){sheet_i}.(exp_list{exp_i}).cell_f_filt;
            cell_f_std(end+1,:) = std(t);
            cell_f_max(end+1,:) = max(t);

            cell_ID_mat(end+1,:) =  All_data2.(varnames{var_i}){sheet_i}.(exp_list{exp_i}).cell_identities;
            var_exp_used{end+1,1} = varnames{var_i};
            var_exp_used{end,2} = exp_list{exp_i};
        end
    end
    
    
    
    % get a list of combined cells
    cell_IDs_comb_t = All_data2.(varnames{var_i}){sheet_i}.(exp_list{exp_i}).cell_IDs_comb;
    cell_IDs_comb_t_unique = [];
    cell_count = 0;
    cells_done = [];
    for cell_i = 1: size(cell_IDs_comb_t,2)
        cells_to_take = cell_IDs_comb_t{cell_i};
        remove_overlapping = find(ismember(cells_to_take,cells_done));
        cells_to_take2 = cells_to_take;
        cells_to_take2(remove_overlapping) = [];
        if isempty(cells_done) || ~any(ismember(cells_done,cells_to_take))
            cell_count = cell_count+1;
            cell_IDs_comb_t_unique{cell_count} = cell_IDs_comb_t{cell_i};
            cells_done = [cells_done,cell_IDs_comb_t{cell_i}];
        elseif ~isempty(cells_to_take2)

            cell_count = cell_count+1;
            cell_IDs_comb_t_unique{cell_count} = cell_IDs_comb_t{cell_i};
            cells_done = [cells_done,cell_IDs_comb_t{cell_i}];
        end
    end
    
    %
    for cell_i = 1:size(cell_ID_mat,2)
        
        cell_list = [cell_ID_mat{:,cell_i}];
        
        [unique_cells] = unique(cell_list,'stable'); % assume that ROIs need to have been common more than once for it to really have been the same cell
        
        c = arrayfun(@(x)length(find(cell_list == x)), unique(cell_list), 'Uniform', false);
        cell_counts = cell2mat(c);
      
        idx_to_del = find(cell_counts<2);
        idx_to_del2 = cellfun(@(c)(cell_list == c),num2cell(unique_cells(find(cell_counts<2))),'UniformOutput',false);
        idx_to_del2 = find(sum(vertcat(idx_to_del2{:}),1) == 1);
        
        unique_cells(idx_to_del) = []; cell_list(idx_to_del2) = [];
        
        cell_IDs_comb(cell_i) = {unique_cells};
    end
    for cell_i = 1:size(classes_mat,1)

        classes_tbl.cell_identities{cell_i} = cell_IDs_comb{cell_i};

    end
    
    
    
    
%     %% Get cell placement info
    clear classes_cmb_cells_tbl
    
    for cell_i = 1:size(cell_IDs_comb_t_unique,2)
%         t = classes_tbl.cell_identities{cell_i};
%         if length(t)>1
        cells_to_take2 = cell_IDs_comb_t_unique{cell_i};    

            sections_t = classes_tbl.sections(cells_to_take2,:);
            maps_t = classes_tbl.maps(cells_to_take2,:);
            X_t = classes_tbl.X(cells_to_take2,:);
            Y_t = classes_tbl.Y(cells_to_take2,:);
           
            exp =  length(maps_t(1,:));
            maps_tt = [];maps_ttt = [];X_tt = [];Y_tt = [];
            for cell_c = 1:length(cells_to_take2)
                X_tt(:,:,cell_c) = X_t{cell_c};
                Y_tt(:,:,cell_c) = Y_t{cell_c};
            end
            for exp_i = 1:exp
                for cell_c = 1:length(cells_to_take2)
                    maps_tt(:,:,cell_c) = maps_t{cell_c,exp_i};
                end
                maps_ttt{exp_i} =  mean(maps_tt,3);
                % figure; imagesc(maps_ttt{1});figure; imagesc(maps_t{1,1}); figure; imagesc(maps_t{2,1});
            end
            sections_t = mean(sections_t);
            maps_t = maps_ttt;
            X_t = mean(X_tt,3);
            Y_t = mean(Y_tt,3);
            % add to the table

            classes_cmb_cells_tbl.cell_IDs(cell_i) = [cell_i];
            if isstruct(classes_cmb_cells_tbl)
                classes_cmb_cells_tbl = struct2table(classes_cmb_cells_tbl);
            end
            classes_cmb_cells_tbl.cell_identities{cell_i} = [cells_to_take2];
            
            % get the new classification and use it if it works
            
%             t = classes_tbl.cell_class(cells_to_take2);
%             unique_classes = unique(t);
%             if length(unique_classes) == 1 % only one class or labelled the same
%                 t = unique_classes;
%             else
%                 
%                 c = arrayfun(@(x)length(find(unique_classes == x)), unique(unique_classes), 'Uniform', false);
%                 class_counts = cell2mat(c);
%                 % check if count is not the same
%                 if length(unique(class_counts)) ~=1 % there is one dominant class
%                     % Take dominant class
%                     dsafsadf
%                 else
%                     % take one with bigger ROI area from the map
%                     map_size = [];
%                     for i = 1:length(cells_to_take2) 
%                         map_size(i) = sum(sum(classes_tbl.maps{cells_to_take2(i),1}));
%                     end
%                     [~,idx] = max(map_size);
%                     t = classes_tbl.cell_class(cells_to_take2(idx));         
%                 end
%             end
%             
% %             classes_cmb_cells_tbl. cell_class{cell_i} = t; old
            classes_cmb_cells_tbl. cell_class{cell_i} = classes_mat2(cell_i,2);
            classes_cmb_cells_tbl. sections_t{cell_i} = sections_t;
            classes_cmb_cells_tbl. maps_t{cell_i} = maps_t;
            classes_cmb_cells_tbl. X_t{cell_i} = X_t;
            classes_cmb_cells_tbl. Y_t{cell_i} = Y_t; % <-this dimension seems to be the scanning dim
            classes_cmb_cells_tbl. Z_t{cell_i} =  sections_t + median(find(mean(Y_t,1)~=0))/length(mean(Y_t,1)) - 1;
            
            classes_cmb_cells_tbl. Zoom{cell_i} =  median([classes_tbl.Zoom{cells_to_take2}]);
            classes_cmb_cells_tbl. Z{cell_i} =  median([classes_tbl.Z{cells_to_take2}]);

    end 
    
    All_data2.classes{which_animal} = classes_tbl;
    All_data2.classes_cmb{which_animal} = classes_cmb_cells_tbl;
%     %% Get maps
    
%     cell_i [groups;classes_cmb_cells_tbl.cell_class_name];
    groupsnos = [groupsnos,classes_cmb_cells_tbl.cell_class{:}];
    coordinate_map_t = [];
        map_all_cells = zeros(size(classes_cmb_cells_tbl.maps_t{1}{1},1),size(classes_cmb_cells_tbl.maps_t{1}{1},2));
        cells_centre = map_all_cells;
        all_x = zeros(size(classes_cmb_cells_tbl.maps_t{1}{1},1),1);
        all_y = zeros(size(classes_cmb_cells_tbl.maps_t{1}{1},1),1);
        for cell_i = 1:height(classes_cmb_cells_tbl)
            groups{end+1} = classes{classes_cmb_cells_tbl. cell_class{cell_i}};
            
            % add all exp maps on top of one another
            mapt = zeros(size(classes_cmb_cells_tbl.maps_t{cell_i}{1},1),size(classes_cmb_cells_tbl.maps_t{cell_i}{1},2),size(classes_cmb_cells_tbl.maps_t{cell_i},2));
%             map_all_cells = map_all_cells + mapt;
            for exp_i = 1:size(classes_cmb_cells_tbl.maps_t{cell_i},2)
                mapt(:,:,exp_i) = classes_cmb_cells_tbl.maps_t{cell_i}{exp_i};
            end
            
             % rotate to make A-P;M-L consistent
             if any(strcmpi(animal_list,sheets{sheet_i})) % ETH experiments
                 mapt = imrotate(mapt,45,'bicubic','crop'); 

                 % make mice face to the left <<<
                 % If stage is moved left  < towards animals P, FoW moves left,                     GRIN flip = A. <     
                 % If stage is moved right > towards animals A, FoW moves right ,                   GRIN flip = P. >  
                 % If stage is moved Up  ^   towards animals L (left hem. implant), FoW moves up ,  GRIN flip = M. ^   
                 % If stage is moved down v towards animals M, FoW moves down ,                     GRIN flip = L. v   
             else
                 mapt = imrotate(mapt,180,'bicubic','crop'); % already

             end
%              % check rotation
%              figure; subplot(1,2,1); imagesc(mapt); subplot(1,2,2); imagesc(imrotate(mapt,45,'crop'));
             % for each cell find center and only take center coordinate
             map_all_cells = map_all_cells+mapt(:,:,1);% for all images will take the first exp only
             
             % find average cell center (median) for all expriments
             tx =[];ty = [];
             for exp_i = 1:size(classes_cmb_cells_tbl.maps_t{cell_i},2)
                [m,iy] = max(mean(mapt(:,:,exp_i),1,'omitnan')); % average of x axis to find max y position
                [m,ix] = max(mean(mapt(:,:,exp_i),2,'omitnan'));
                
                ix = round(median(ix)); % maybe more than 1!
                iy = round(median(iy));
                tx(exp_i) = ix;
                ty(exp_i) = iy;
             end
             ix = round(median(tx)); % find center
             iy = round(median(ty));
             iz = classes_cmb_cells_tbl.Z_t{cell_i}/5;
             
            cells_centre(ix,iy) = cells_centre(ix,iy)+1;
            
            %%% adjust coordinates to um
            pxl_no = size(mapt,1);
            um_plane = 600; 
            zoom_t = classes_cmb_cells_tbl.Zoom{cell_i};
            z = (classes_cmb_cells_tbl.Z{cell_i}/2); % dividing by 2 for pointspread correction; dividing by plane number
            
            um_plane = (um_plane * um_plane)/(um_plane*zoom_t);
            ratio_plane = um_plane/pxl_no;
            
   
            coordinate_map_t(cell_i,1) = ix*ratio_plane;
            coordinate_map_t(cell_i,2) = iy*ratio_plane;
            coordinate_map_t(cell_i,3) = iz*z;%section_i;
            
%              [m,ix] = max(median(mean(mapt,1,'omitnan'),3,'omitnan'));
%              [m,iy] = max(median(mean(mapt,2,'omitnan'),3,'omitnan'));
%              cells_centre(ix,iy) = cells_centre(ix,iy)+1;
             %             all_x = all_x + all_cells_in_section.X{cell_i}';
             %             all_y = all_y + all_cells_in_section.Y{cell_i};
%               figure; imshowpair(mapt(:,:,3:5),cells_centre);
%              fsadf
        end

%         figure; imagesc(mapt(:,:,1:3)+cells_centre); %hold on; imagesc(cells_centre);
        
        
        
        %%% center cells on the axes system
        
        center_x_y_z(1) = max(coordinate_map_t(:,1))/2;
        center_x_y_z(2) = max(coordinate_map_t(:,2))/2;
        center_x_y_z(3) = max(coordinate_map_t(:,3))/2;
%         center_x_y_z(3) = 0;
%         center_x_y_z = [0,0,0];

        coordinate_map_t(:,1) =  coordinate_map_t(:,1) - center_x_y_z(1);
        coordinate_map_t(:,2) =  coordinate_map_t(:,2) - center_x_y_z(2);
        coordinate_map_t(:,3) =  coordinate_map_t(:,3) - center_x_y_z(3);
        
        if size(coordinate_map_t,1) ~= height(classes_cmb_cells_tbl)
            safsafddsf
        end
        coordinate_map = [coordinate_map;coordinate_map_t];
        coordinate_map_animal{sheet_i} = coordinate_map_t;
%         
        all_cells_map{sheet_i} = map_all_cells;
%         groups_animal{sheet_i,section_i,1} = classes_cmb_cells_tbl.cell_class_name;
        groups_animal{sheet_i,2} = [classes_cmb_cells_tbl.cell_class{:}];
           
%          s = scatterhist(coordinate_map(:,1),coordinate_map(:,2),'Kernel','on','Group',groups);
%          s.YLimits
%          s.YLimits = [1 256];
%          s.XLimits = [1 256];
         
%          xlim([1,256]);ylim([1,256]);
%           axis square
%         set(gca,'xtick',[],'ytick',[]);
%     end
end
%% Look at running type
animal_list = {'077','082'}; % OK Running speed
animal_list2 = {'126b','152f','G35','G36','m143a'}; % too fast - is x10?

% For each cell, each exp, runn xcorr vs running -> Get exp number of
% corrs, per cell. Mean for corr score for each cell.
% 1) Look if any glu classes are more represented in corr scores
% 2)Check for anatomical distribution of xcor scores

for animal_i = 1:height(All_data2)
    cells_animal = height(All_data2.classes_cmb{animal_i});
    
    for var_i = 1:length(varnames)
        exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
        for exp_i = 1:length(exp_list)
            
            fps = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).fps/6;
            cell_sig_t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common(:,1);
            run_t =  All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).abs_run_vec;
            dif_sampling = length(run_t)/length(cell_sig_t);
            run_t = resample(run_t,1*10000,round(dif_sampling*10000));
            
            if length(cell_sig_t) ~= length(run_t)
                dif_sampling = length(run_t)/length(cell_sig_t);
                run_t = resample(run_t,1*10000,round(dif_sampling*10000));
            end
            for cell_i = 1:cells_animal % check all experiments for each cell and get their scores
                cell_sig_t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common(:,cell_i);
                [Spearman_t, Spearman_p] = corr(cell_sig_t, run_t, 'type', 'Spearman');
                [Pearson_t, Pearson_p]  = corr(cell_sig_t, run_t, 'type', 'Pearson');
                xcorr(cell_sig_t, run_t);
                [c,l]  = xcorr(cell_sig_t,run_t,'coeff'); % positive means first vector is laging behind. Negative means 2nd vector is behind in time.
                l = l/fps;
%                 plot(l,c);hold on; xlabel('Lag time (s)');ylabel('Correlation score'); xline(0);
                
                max_c = max(c);
                lag = l(min(abs(find(c == max(c)))))/fps;
                
                
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Spearman_t')){cell_i,exp_i} = Spearman_t;
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Spearman_p')){cell_i,exp_i} = Spearman_p;
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Pearson_t')){cell_i,exp_i} = Pearson_t;
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Pearson_p')){cell_i,exp_i} = Pearson_p;
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Xcor_maxC')){cell_i,exp_i} = max_c;
                All_data2.classes_cmb{animal_i}.(strcat(varnames{var_i},'_Xcor_Lag')){cell_i,exp_i} = lag; % seconds
            end
        end
    end
end
%% Check for each class type, plot mean + sem of a given corr score
Corr_score_mat = [];
close all
% hold on
for animal_i = 1:height(All_data2) % [1:3,6]%
   
    t = All_data2.classes_cmb{animal_i}.Sal_Pearson_t;
    p_idx = All_data2.classes_cmb{animal_i}.Sal_Pearson_p;
    
%     t = All_data2.classes_cmb{animal_i}.Sal_Spearman_t;
%     p_idx = All_data2.classes_cmb{animal_i}.Sal_Spearman_p;
    p_idx = cell2mat(p_idx);
    p_idx = mean(p_idx,2);
    t = cell2mat(t);
   t = mean(t,2);
   p_idx = find(p_idx<1000); % all
%     p_idx = find((p_idx < 0.05 & t > 0.1) | (p_idx < 0.05 & t< -0.1));
%     p_idx = find((p_idx > 0.05 & t < 0.1) | (p_idx > 0.05 & t> -0.1));
%    t = All_data2.classes_cmb{animal_i}.Glu_Pearson;
%  t = All_data2.classes_cmb{animal_i}.Sal_Xcor_maxC;
%  t = All_data2.classes_cmb{animal_i}.Sal_Xcor_Lag;
   
   t = t(p_idx);
   types = All_data2.classes_cmb{animal_i}.cell_class;
   height_t = size(Corr_score_mat,1);
   Corr_score_mat(height_t+1:height_t+size(t,1),1) = cell2mat(types(p_idx));
   Corr_score_mat(height_t+1:height_t+size(t,1),2) = t;
   Corr_score_mat(height_t+1:height_t+size(t,1),3) = animal_i;
end


types = unique(Corr_score_mat(:,1));
tmat = nan(10000,length(types));
count_cell_types = [];
for type_i = 1:length(types)
    idx = find(Corr_score_mat(:,1) == types(type_i));
    tmat(1:length(idx),type_i) = Corr_score_mat(idx,2);
%     m_t = mean(Corr_score_mat(:,idx),'omitnan');
%     sem_t = nansem_nansuite(Corr_score_mat(:,idx));
    count_cell_types(type_i) = sum(Corr_score_mat(:,1)==type_i);
end
% %%
notBoxPlot_mine(tmat,'style','errorbar.');
xtick = [1:5];
xticklabels(char(classes));
% [P,ANOVATAB,STATS] =anova1(tmat);
figure;
[p,tbl,stats] = anova1(tmat,[] ,'off'); % ,'off'
stats.gnames = classes';
[c,m,h,nms] = multcompare(stats); %
p_report =   round(c(:,end),3)';


figure;
xPercent = count_cell_types / sum(count_cell_types) * 100;
newLabels = [];
for i=1:length(count_cell_types)
%     newLabels = [newLabels {sprintf('%s %i/%i (%.1f%%)', labels{i}, cnt_unique(i),sum(cnt_unique),xPercent(i))}]; 
    newLabels = [newLabels {sprintf('%s %.1f%%', classes{i}, xPercent(i))}];   
end

h1 = pie(count_cell_types,newLabels);
ax = gca(); 
ax.Colormap = colorlist(1:length(count_cell_types),:); 


% Create legend
% labels = {'G','iG','dG','idG','None'};
labels = classes;
legend(labels);
h1 = legend;
h1.FontSize = 6;
h1.Location = "southeastoutside";
% title('Positive correlation');
title(strjoin({'Cells grouped by glucose response n =', num2str(total_num_cells)}));
%% Spit into correlated and non-correlated cells
Corr_score_mat = []; Pos_Corr_score_mat = []; Neg_Corr_score_mat = []; No_Corr_score_mat = [];
close all
% hold on
    cut_off_critearea = [0.05,0.01]; % p value t
for animal_i = 1:height(All_data2) % [1:3,6]%
    
    t_keep = All_data2.classes_cmb{animal_i}.Sal_Spearman_t; %Sal_Spearman_t Sal_Pearson_t
    p_keep = All_data2.classes_cmb{animal_i}.Sal_Spearman_p; % Sal_Spearman_p Sal_Pearson_p
    
    % Postive corr
    p_temp = cut_off_critearea(1);
    m = size(t_keep,1); % null hypotheses number
    p_temp = 1 - (1-p_temp)^(1/m);
    
    p_idx = cell2mat(p_keep);
    p_idx = mean(p_idx,2);
    t = cell2mat(t_keep);
    t = mean(t,2);
    p_idx = find((p_idx < p_temp & t > cut_off_critearea(2)));
    p_idx_keep = p_idx;

    
    t = t(p_idx);
    types = All_data2.classes_cmb{animal_i}.cell_class;
    height_t = size(Pos_Corr_score_mat,1);
    Pos_Corr_score_mat(height_t+1:height_t+size(t,1),1) = cell2mat(types(p_idx));
    Pos_Corr_score_mat(height_t+1:height_t+size(t,1),2) = t;
    Pos_Corr_score_mat(height_t+1:height_t+size(t,1),3) = animal_i;
    
    % Negative corr
    p_temp = cut_off_critearea(1);
    m = size(t_keep,1); % null hypotheses number
    p_temp = 1 - (1-p_temp)^(1/m);
    
    p_idx = cell2mat(p_keep);
    p_idx = mean(p_idx,2);
    t = cell2mat(t_keep);
    t = mean(t,2);
    p_idx = find((p_idx < p_temp & t < -cut_off_critearea(2)));
    p_idx_keep = [p_idx_keep;p_idx];
    
    t = t(p_idx);
    types = All_data2.classes_cmb{animal_i}.cell_class;
    height_t = size(Neg_Corr_score_mat,1);
    Neg_Corr_score_mat(height_t+1:height_t+size(t,1),1) = cell2mat(types(p_idx));
    Neg_Corr_score_mat(height_t+1:height_t+size(t,1),2) = t;
    Neg_Corr_score_mat(height_t+1:height_t+size(t,1),3) = animal_i;
    
    
    % No corr
    p_temp = cut_off_critearea(1);
    m = size(t_keep,1); % null hypotheses number
    p_temp = 1 - (1-p_temp)^(1/m);
    
    p_idx = cell2mat(p_keep);
    p_idx = mean(p_idx,2);
    t = cell2mat(t_keep);
    t = mean(t,2);
%     p_idx = find((p_idx > p_temp |  p_idx < p_temp & t < cut_off_critearea(2)& t > -cut_off_critearea(2)) ) ;
    

    p_idx = [1:length(t)];
    p_idx(p_idx_keep) = [];
    t = t(p_idx);
    types = All_data2.classes_cmb{animal_i}.cell_class;
    height_t = size(No_Corr_score_mat,1);
    No_Corr_score_mat(height_t+1:height_t+size(t,1),1) = cell2mat(types(p_idx));
    No_Corr_score_mat(height_t+1:height_t+size(t,1),2) = t;
    No_Corr_score_mat(height_t+1:height_t+size(t,1),3) = animal_i;
    
end

total_num_cells = height(No_Corr_score_mat)+height(Neg_Corr_score_mat)+height(Pos_Corr_score_mat);
Corr_score_mat = No_Corr_score_mat;
types = unique(Corr_score_mat(:,1));
tmat = nan(10000,length(types));
for type_i = 1:length(types)
    idx = find(Corr_score_mat(:,1) == types(type_i));
    tmat(1:length(idx),type_i) = Corr_score_mat(idx,2);
%     m_t = mean(Corr_score_mat(:,idx),'omitnan');
%     sem_t = nansem_nansuite(Corr_score_mat(:,idx));    
end
% %%
notBoxPlot_mine(tmat,'style','errorbar.');
xtick = [1:5];
xticklabels(char(classes));
% [P,ANOVATAB,STATS] =anova1(tmat);

[p,tbl,stats] = anova1(tmat,[] ,'off'); % ,'off'
stats.gnames = classes';%{'G','iG','dG','idG','None'}';
[c,m,h,nms] = multcompare(stats); %
p_report =   round(c(:,end),3)';

% % pie charts for each cell type
% figure('units','normalized','outerposition',[0 0 0.5 1]);


% pie charts for each condition
figure('units','normalized','outerposition',[0 0 0.5 1]);
set(gcf,'color','w');
subplot(2,2,1);
a = Pos_Corr_score_mat(:,1)';

[cnt_unique, unique_a] = hist(a,unique(a),'off');
explode = [0 1 0 1 0];
explode = [1 1 1 1 1];
explode = [0 0 0 0 0];
% h1 = pie(cnt_unique,explode);

labels = classes;%{'G','iG','dG','idG','None'};
xPercent = cnt_unique / sum(cnt_unique) * 100;
newLabels = [];
for i=1:length(cnt_unique)
%     newLabels = [newLabels {sprintf('%s %i/%i (%.1f%%)', labels{i}, cnt_unique(i),sum(cnt_unique),xPercent(i))}]; 
    newLabels = [newLabels {sprintf('%s %.1f%%', labels{i}, xPercent(i))}];   
end

h1 = pie(cnt_unique,explode,newLabels);

ax = gca(); 
ax.Colormap = colorlist(1:length(unique_a),:); 


% Create legend
% labels = {'G','iG','dG','idG','None'};
labels = labels(unique_a);
legend(labels);
h1 = legend;
h1.FontSize = 6;
h1.Location = "southeastoutside";
% title('Positive correlation');
title(strjoin({'Positive correlation cells', num2str(sum(cnt_unique)),'/',num2str(total_num_cells)}));

subplot(2,2,2);
a = Neg_Corr_score_mat(:,1)';
[cnt_unique, unique_a] = hist(a,unique(a),'off');
labels = classes;%{'G','iG','dG','idG','None'};
xPercent = cnt_unique / sum(cnt_unique) * 100;
newLabels = [];
for i=1:length(cnt_unique)
%     newLabels = [newLabels {sprintf('%s %i/%i (%.1f%%)', labels{i}, cnt_unique(i),sum(cnt_unique),xPercent(i))}];   
    newLabels = [newLabels {sprintf('%s %.1f%%', labels{i}, xPercent(i))}];
end

h1 = pie(cnt_unique,newLabels);
ax = gca(); 
ax.Colormap = colorlist(1:length(unique_a),:); 

% Create legend
% labels = {'G','iG','dG','idG','None'};
% labels = labels(unique_a);
legend(labels);
h1 = legend;
h1.FontSize = 6;
h1.Location = "southeastoutside";
% title('Negative correlation');
title(strjoin({'Negative correlation cells', num2str(sum(cnt_unique)),'/',num2str(total_num_cells)}));

subplot(2,2,3);
a = No_Corr_score_mat(:,1)';
[cnt_unique, unique_a] = hist(a,unique(a),'off');

labels = classes;%{'G','iG','dG','idG','None'};
xPercent = cnt_unique / sum(cnt_unique) * 100;
newLabels = [];
for i=1:length(cnt_unique)
%     newLabels = [newLabels {sprintf('%s %i/%i (%.1f%%)', labels{i}, cnt_unique(i),sum(cnt_unique),xPercent(i))}];  
    newLabels = [newLabels {sprintf('%s %.1f%%', labels{i}, xPercent(i))}];
end

h1 = pie(cnt_unique,newLabels);
ax = gca(); 
ax.Colormap = colorlist(1:length(unique_a),:); 

% Create legend

labels = labels(unique_a);
legend(labels);
h1 = legend;
h1.FontSize = 6;
h1.Location = "southeastoutside";
title(strjoin({'No correlation cells', num2str(sum(cnt_unique)),'/',num2str(total_num_cells)}));
suptitle(strjoin({'Split by correlations with running',newline,'p <',num2str(cut_off_critearea(1)), 'rho >',num2str(cut_off_critearea(2))}));

% pie charts for each state
state_mat = [];
state_mat(1:size(Pos_Corr_score_mat,1),:) = Pos_Corr_score_mat;
state_mat(1:size(Pos_Corr_score_mat,1),4) = 1;
state_mat(size(state_mat,1)+1:size(state_mat,1)+size(Neg_Corr_score_mat,1),1:3) = Neg_Corr_score_mat;
state_mat(end-size(Neg_Corr_score_mat,1)+1:end,4) = 2;
state_mat(size(state_mat,1)+1:size(state_mat,1)+size(No_Corr_score_mat,1),1:3) = No_Corr_score_mat;
state_mat(end-size(No_Corr_score_mat,1)+1:end,4) = 3;

figure('units','normalized','outerposition',[0.5 0 0.5 1]);
set(gcf,'color','w');

labels = classes;%{'G','iG','dG','idG','None'};

% Create legend
% labels2 = {'Pos cor','Neg cor','No cor'};

% labels = labels(unique_a);


for state_i = 1:length(labels);
    if state_i > 2
        subplot(2,3,state_i+1);
    else
        subplot(2,3,state_i);
    end
    a = [];
    a = state_mat(find(state_mat(:,1) == state_i),4);
    [cnt_unique, unique_a] = hist(a,unique(a),'off');
    
    labels2 = {'Pos cor','Neg cor','No cor'};
    xPercent = cnt_unique / sum(cnt_unique) * 100;
    newLabels = [];
    for i=1:length(cnt_unique)
%         newLabels = [newLabels {sprintf('%s %i/%i (%.1f%%)', labels2{i}, cnt_unique(i),sum(cnt_unique),xPercent(i))}];
        newLabels = [newLabels {sprintf('%i/%i (%.1f%%)', cnt_unique(i),sum(cnt_unique),xPercent(i))}];
    end

    h1 = pie(cnt_unique,newLabels);
%     h1 = pie(cnt_unique);
    ax = gca();
    ax.Colormap = colorlist(1:length(unique_a),:);
    
    
%     legend(labels2);
%     h1 = legend;
%     h1.FontSize = 6;
%     h1.Location = "southeastoutside";
    title(labels(state_i));
end


% suptitle('Split by correlations with glucose response');
suptitle(strjoin({'Split by correlations with glucose response',newline,'p <',num2str(cut_off_critearea(1)), 'rho >',num2str(cut_off_critearea(2))}));
%% Combine glucose inhibited cells
comb_idxes = find(strcmpi(classes,'iG')| strcmpi(classes,'idG'));

figure;
  a = [];

  
  
    a = state_mat(find(state_mat(:,1) == comb_idxes(1) | state_mat(:,1) == comb_idxes(2)),4);
    
    [cnt_unique, unique_a] = hist(a,unique(a),'off');
    h1 = pie(cnt_unique);
    ax = gca();
    ax.Colormap = colorlist(1:length(unique_a),:);
    
    
    legend(labels2);
    h1 = legend;
    h1.FontSize = 6;
    h1.Location = "southeastoutside";
    

title('Glucose inhibited cells');
%%
close all

colorlist = [0,0,0;100,143,255;120,94,240;220,38,127;254,97,0;255,176,0]/255;
% colorlist = [0,13,255;40,120,255;255,30,30;255,100,0;70,62,70;0,0,0]/255;
% sort by classes
[~,sort_classes] = sort(groupsnos);

coordinate_map = coordinate_map(sort_classes,:);
groups = groups(sort_classes);
if size(groups,2)> size(groups,1)
    groups = groups';
end
groupsnos = groupsnos(sort_classes);

figure;

 h1 = scatterhist(coordinate_map(:,1),coordinate_map(:,2),'Kernel','on','Group',groups,'Color',colorlist(1:length(unique(groups)),:),'Marker','x','LineStyle','-');
 %   set(gca,'xtick',[],'ytick',[]);
 xlabel('A->P');
 ylabel('M->L');
 title('All animals');
 
 
%  hold on;
%  clr = get(h1(1),'colororder');
%  boxplot(h1(2),coordinate_map(:,1)',groups,'orientation','horizontal',...
%      'color',clr);
%  boxplot(h1(3),coordinate_map(:,2)',groups,'orientation','horizontal',...
%      'color',clr);
 
%   %% X vs Z
  figure;
 h2 =scatterhist(coordinate_map(:,1),coordinate_map(:,3),'Kernel','on','Group',groups,'Color',colorlist(1:length(unique(groups)),:),'Marker','x','LineStyle','-');
%   set(gca,'xtick',[],'ytick',[]);
   xlabel('A->P');
  ylabel('V->D');
  title('All animals');
  
%  hold on;
%  clr = get(h2(1),'colororder');
%  boxplot(h2(2),coordinate_map(:,1)',groups,'orientation','horizontal',...
%      'color',clr);
%  boxplot(h2(3),coordinate_map(:,3)',groups,'orientation','horizontal',...
%      'color',clr);
%  
%     %% X vs Z
  figure;
 h3 = scatterhist(coordinate_map(:,2),coordinate_map(:,3),'Kernel','on','Group',groups,'Color',colorlist(1:length(unique(groups)),:),'Marker','x','LineStyle','-');
%   set(gca,'xtick',[],'ytick',[]);
   xlabel('M->L');
  ylabel('V->D');
  title('All animals');
  
%  hold on;
%  clr = get(h3(1),'colororder');
%  boxplot(h3(2),coordinate_map(:,2)',groups,'orientation','horizontal',...
%      'color',clr);
%  boxplot(h3(3),coordinate_map(:,3)',groups,'orientation','horizontal',...
%      'color',clr);
 
  %%
  figure;
  
  [c2,ia2,ic2]  = unique(groups);
  
  % color cells accordingly
  ic2_color = zeros(length(ic2),3);
  for i = 1:length(ic2)
      ic2_color(i,:) = colorlist(ic2(i,:));
  end

 hh = scatter3(coordinate_map(:,1),coordinate_map(:,2),coordinate_map(:,3),ones(length(groupsnos),1)*50,groupsnos,'Marker','x');
 xlabel('A->P');
 ylabel('M->L');
 zlabel('V->D');
 
 figure; 
 scatter3([1:5],[1:5],[1:5],ones(5,1)*50,[1:5],'Marker','x'); title('colormap');
%% For each class find distributions over V->D axis
% close all
% [~,sort_classes] = sort(groupsnos);
% coordinate_map = coordinate_map(sort_classes,:);
% groups = groups(sort_classes);
% 
% groupsnos = groupsnos(sort_classes);
axes = {'A->P','M->L','V->D'};
for axis_i = 1:3
    groups_list = [1:length(unique(groupsnos))];
    lgnd = unique(groups,'stable');
    coordinates_of_group = nan(length(groups),length(unique(groupsnos)));
    for group_i = groups_list
        idx_of_group = find(groupsnos == groups_list(group_i));
        coordinates_of_group(idx_of_group,group_i) = coordinate_map(idx_of_group,axis_i);
    end
    figure;
    notBoxPlot_mine(coordinates_of_group,'style','errorbar.');
    ylabel(strjoin({axes{axis_i},'(um)'}));
    % XTickLabel()
    % anova1(coordinates_of_group)
    
    [p,tbl,stats] = anova1(coordinates_of_group,[] ,'off'); % ,'off'
    [c,m,h,nms] = multcompare(stats,[] ,'off'); %
    p_report =   round(c(:,end),3)';
    title(strjoin({'Axis =',axes{axis_i},newline,'ANOVA p =',num2str(round(p,5))}));   xticklabels(lgnd);
    % figure;[c,m,h,nms] = multcompare(stats,[] );yticklabels(flip(lgnd));%
end
 
% %% Plot per animal
% % close all
% % plot per animal to check
% for animal_i = 1:size(coordinate_map_animal,2)
%     animal_coordinate_map = [];
%     group_no_t = [];
%     group_name_t = {};
% %     for section_i = 1:size(coordinate_map_animal,2)
%         % combine sections per animal
%         animal_coordinate_map = [animal_coordinate_map;coordinate_map_animal{animal_i}];
%         
%         group_no_t = groups_animal{animal_i,2};
%         for cell_i = 1:length(group_no_t)
%             group_name_t{cell_i} = classes{group_no_t(cell_i)};
%         end
% %     end
%     [~,sort_classes] = sort(group_no_t);
%     animal_coordinate_map = animal_coordinate_map(sort_classes,:);
%     group_name_t = group_name_t(sort_classes)';
% 
%     
% %     figure;
% %     subplot(3,1,1);
% %     scatterhist(animal_coordinate_map(:,1),animal_coordinate_map(:,2),'Kernel','on','Group',group_name_t,'Color',colorlist(1:length(unique(group_name_t)),:),'Marker','x','LineStyle','-');
% % %     set(gca,'xtick',[],'ytick',[]);
% %     xlabel('A->P');
% %     ylabel('M->L');
% %     title(sheets{animal_i});
%     %   %% X vs Z
% %     figure;
% %     scatterhist(animal_coordinate_map(:,1),animal_coordinate_map(:,3),'Kernel','on','Group',group_name_t,'Color',colorlist(1:length(unique(group_name_t)),:),'Marker','x','LineStyle','-');
% %     set(gca,'xtick',[],'ytick',[]);
% %     xlabel('A->P');
% %     ylabel('V->D');
% %     title(sheets{animal_i});
% %     
% %     %     %% X vs Z
%     figure;
%     scatterhist(animal_coordinate_map(:,2),animal_coordinate_map(:,3),'Kernel','on','Group',group_name_t,'Color',colorlist(1:length(unique(group_name_t)),:),'Marker','x','LineStyle','-');
% %     set(gca,'xtick',[],'ytick',[]);
%     xlabel('M->L');
%     ylabel('V->D');
%     title(sheets{animal_i});
% %     
% end


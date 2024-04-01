% Load pre-processed 2P sal vs glu data into a single struct
%%%%%% Analysis workflow:
% cell_f_extraction_2P_doublecheckplanes_24012022a.m
% ExportPreprocesseddata_2PGlusensstim08022022.m
% Combine_2Pdata0802.m
% Analysis_combined_data_2P.m
% Locate2Pcells.m

clear all
clc
params.colour_list = {[0 0.4470 0.7410],[ 0.8500 0.3250 0.0980],[ 0.9290 0.6940 0.1250],[ 0.4940 0.1840 0.5560],[ 0.4660 0.6740 0.1880],[ 0.3010 0.7450 0.9330],[ 0.6350 0.0780 0.1840],[0,0,0],[0.5,0,0],[0,0.5,0],[0,0,0.5]};
    
fid = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp';
animal_list = {'077','082'};
animal_list2 = {'126b','152f','G35','G36','m143a'};
animal_list3 = {animal_list{:},animal_list2{:}};
Available_animal_list = getsubfolder_names_fun(fid);


x = cellfun(@(c)strcmpi(Available_animal_list,c),animal_list3,'UniformOutput',false);
x = vertcat(x{:});
x = sum(x,1);

idx_to_del = find(x == 0);
Available_animal_list(:,idx_to_del) = [];
%%
All_data = [];
for animal_i = 1:length(Available_animal_list)
    exp_list = getsubfolder_names_fun(fullfile(fid,Available_animal_list{animal_i}));
     for exp_i = 1:length(exp_list)
        fid_exp = fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i});
        if exist(fullfile(fid_exp,'preprocessed_data.mat')) == 2 %&& 1==2
            disp(strcat('Adding: ',Available_animal_list{animal_i},'-',exp_list{exp_i}));
            load(fullfile(fid,Available_animal_list{animal_i},exp_list{exp_i},'preprocessed_data.mat'));
            All_data(animal_i).ID = Available_animal_list{animal_i};
            if contains(exp_list{exp_i},'sal')
                 All_data(animal_i).Sal.(strcat('exp_',strrep(exp_list{exp_i},'-','_'))) = preprocessed_data;
            else
                All_data(animal_i).Glu.(strcat('exp_',strrep(exp_list{exp_i},'-','_'))) = preprocessed_data;
            end           
        end
     end
end

%% plot mean of all cells sal vs glucose
if isstruct(All_data)
    All_data = struct2table(All_data);
end
varnames = fieldnames(All_data);
varnames = varnames(2:3);

%% Here combine cells that are the same between layers TODO!!
expand_pxls = 3;
overlap_thrs = 0.05; % part of 1

% All_data_processed = [];
% All_data_processed.ID = All_data.ID;
All_data2 = All_data;
for animal_i = 1:height(All_data)
    cell_ID_mat = [];
    plane_mat = [];
    ROI_mat = [];
    var_exp_used = {};
    for var_i = 1:length(varnames)
        
        exp_list = fieldnames(All_data.(varnames{var_i}){animal_i});
        for exp_i = 1:length(exp_list)
            disp(exp_list{exp_i});
            data = All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i});
            try
                cell_f_t = [All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f];
                np_f_t = [All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}).np_f];
            catch
                cell_f_t = [];  np_f_t = []; sizes = [];
                for plane_i = 1:size(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}),2)
                    sizes(plane_i) = length(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_f);
                end
                for plane_i = 1:size(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}),2)
                    All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_f = All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_f(1:min(sizes),:);
                    All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).np_f = All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).np_f(1:min(sizes),:);
                    
                end
                cell_f_t = [All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f];
                np_f_t = [All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}).np_f];
            end
            plane_list = []; ROI_mat_t = {};
            for plane_i = 1:length(data)
                ROI_mat_t = {ROI_mat_t{:},All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).ROIsmat{:}};
                plane_list = [plane_list,ones(1,length(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID))*plane_i];
            end
            
            cell_identities = {};
            % for each cell compare to other cells
            for cell_i = 1:size(cell_f_t,2)
                cell_identities{cell_i} = cell_i;
                for cell_vs = 1:size(cell_f_t,2)
                    if cell_i ~= cell_vs % if not the same cell
                        ROI1 = ROI_mat_t{cell_i};
                        ROI2 = ROI_mat_t{cell_vs};
                        
                        if plane_list(cell_i) == plane_list(cell_vs)% the cell is on the same plane expand ROIs by expand_pxls
                            ROI1_exp = movmean(ROI1,expand_pxls,1);
                            ROI1_exp = movmean(ROI1_exp,expand_pxls,2);
                            ROI1_exp = ceil(ROI1_exp);
                            
                            ROI2_exp = movmean(ROI2,expand_pxls,1);
                            ROI2_exp = movmean(ROI2_exp,expand_pxls,2);
                            ROI2_exp = ceil(ROI2_exp);
                            
                            ROI1 = ROI1_exp;
                            ROI2 = ROI2_exp;
%                             t =imfuse(repmat(ROI1,1,1,3),repmat(ROI1_exp,1,1,3));
%                                imagesc(t)
                        end
                        ROI1v = reshape(ROI1,1,[]);
                        ROI2v = reshape(ROI2,1,[]);
                        
                        idx1 = find(ROI1v == 1);
                        idx2 = find(ROI2v == 1);
                        Overlap_percentages = [length(find(ROI2v(idx1)==1))/length(idx1),length(find(ROI1v(idx2)==1))/length(idx2)];
                     
                        cell_to_take_f = smoothdata(zscore(cell_f_t(:,cell_i)),'movmedian',10);
                        cell_to_compare_to_f = smoothdata(zscore(cell_f_t(:,cell_vs)),'movmedian',10);
                        recorded_framelen = min([length(cell_to_take_f),length(cell_to_compare_to_f)]);
                        cell_to_take_f = cell_to_take_f(1:recorded_framelen);cell_to_compare_to_f = cell_to_compare_to_f(1:recorded_framelen);
                        [C,lags] = xcorr(cell_to_take_f,cell_to_compare_to_f,'coeff');
                        [max_c,lag_at_c] = max(C);
                        lag_at_c = lags(lag_at_c);
                        
                        if any(Overlap_percentages>overlap_thrs) && abs(lag_at_c)<10 && max_c >0.9
                                                     figure;
                                                     subplot(2,1,1);
                                                     t =imfuse(repmat(ROI1,1,1,3),repmat(ROI2,1,1,3));
                                                     imagesc(t)
                                                     subplot(2,1,2);
                                                     plot(cell_to_take_f); hold on; plot(cell_to_compare_to_f);
                                                     title(strjoin({'C =',num2str(max_c),'Lag =',num2str(lag_at_c),'Overlap =',num2str(max(Overlap_percentages))}));
                                                    suptitle(strjoin({'Plane',num2str(plane_i),', cell,',num2str(cell_i),'&',num2str(cell_vs)}));
                             cell_identities{cell_i} = [ cell_identities{cell_i},cell_vs];
% cell_similar_to(end+1) = cell_vs;
                        end
                    end
                end
            end 
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).cell_f = cell_f_t;
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).np_f = np_f_t;
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).plane_list = plane_list;
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).cell_identities = cell_identities;
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).ROI_mat = ROI_mat_t;
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(2:5) = [];
        end
        
    end
end
save(fullfile(fid,'All_data.mat'),'All_data2', '-v7.3');
     
%% Filter signal subtracting noise (faster than 5Hz);
for animal_i = 1:height(All_data2)
    for var_i = 1:length(varnames)
        exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
        plane_list = All_data2.(varnames{var_i}){animal_i}.(exp_list{1}).plane_list;
        for exp_i = 1:length(exp_list)
            for cell_i = 1:size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f,2)
                t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f(:,cell_i);
                tt = t;
                padding_size = 500;
                tt = [repmat(t(1,:),padding_size, 1); t ; repmat(t(end,:),padding_size, 1)]; % pad for correct resampling
                tt = highpass(tt,0.2,5.2,'Steepness',0.8); % resample to the given frequency
                tt = tt(padding_size+1:end-padding_size); % remove padding
%                 plot(t-tt);
                All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt(:,cell_i) = t-tt;
            end
        end
    end
end
save(fullfile(fid,'All_data.mat'),'All_data2', '-v7.3');
%% Combine ROIs for cells that have been identified as identical:
% 1) if ROI was identified in all sessions = mean to get cell
% 2) If discrepancies between sessions occur, take ROI with highest
% fluroescence STD
% 2B) Mean all ROIs for a given cell
for animal_i = 1:height(All_data2)
    cell_ID_mat = {};
    cell_IDs_comb = {};
    cell_f_std = [];
    cell_f_max = [];
    var_exp_used = {};
    for var_i = 1:length(varnames)
        
        exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
        plane_list = All_data2.(varnames{var_i}){animal_i}.(exp_list{1}).plane_list;
        count = 1;
        for exp_i = 1:length(exp_list)

            t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt;
            cell_f_std(end+1,:) = std(t);
            cell_f_max(end+1,:) = max(t);

            cell_ID_mat(end+1,:) =  All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_identities;
            var_exp_used{end+1,1} = varnames{var_i};
            var_exp_used{end,2} = exp_list{exp_i};
        end
    end
%     cell_f_stdmax = cell_f_std./cell_f_max;

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
    
    for var_i = 1:length(varnames)
        
        exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
        
        plane_list = All_data2.(varnames{var_i}){animal_i}.(exp_list{1}).plane_list;
        count = 1;
        for exp_i = 1:length(exp_list)
            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_IDs_comb = cell_IDs_comb;
            cell_IDs_comb_t = cell_IDs_comb;
            
            
            sorted_stimframes = round(sort(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).stim_frames_to_use)/6);
            last_stim_frame = sorted_stimframes(length(sorted_stimframes)/2-1);
            start_frame = 1*60* round(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).fps / 6);
            post_IP_frame = sorted_stimframes(length(sorted_stimframes)/2+1);

            All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common = [];
            cell_count = 0;
            cells_done = [];
            for cell_i = 1: size(cell_IDs_comb_t,2)
                cells_to_take = cell_IDs_comb{cell_i};
                remove_overlapping = find(ismember(cells_to_take,cells_done));
                cells_to_take2 = cells_to_take;
                cells_to_take2(remove_overlapping) = [];
                if isempty(cells_done) || ~any(ismember(cells_done,cells_to_take))
                    t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt(:,cell_IDs_comb{cell_i});
                    tt = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f(:,cell_IDs_comb{cell_i});
                    for cell_ii = 1:size(t,2)
                        m = mean(t(start_frame:last_stim_frame,cell_ii),'omitnan');
                        std_t = std(t(start_frame:last_stim_frame,cell_ii),'omitnan');
                        t(:,cell_ii) = (t(:,cell_ii) - m)/std_t;
                        
                        m = mean(tt(start_frame:last_stim_frame,cell_ii),'omitnan');
                        std_t = std(tt(start_frame:last_stim_frame,cell_ii),'omitnan');
                        tt(:,cell_ii) = (tt(:,cell_ii) - m)/std_t;
                    end
                    %                 tt = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f(:,cell_IDs_comb{cell_i});
                    cell_count = cell_count+1;
                    All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common(:,cell_count) = mean(t,2,'omitnan');
                    All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common(:,cell_count) = mean(tt,2,'omitnan');
                    cells_done = [cells_done,cell_IDs_comb{cell_i}];
                elseif ~isempty(cells_to_take2)
                    t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt(:,cells_to_take2);
                    tt = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f(:,cells_to_take2);
                    for cell_ii = 1:size(t,2)
                        m = mean(t(start_frame:last_stim_frame,cell_ii),'omitnan');
                        std_t = std(t(start_frame:last_stim_frame,cell_ii),'omitnan');
                        t(:,cell_ii) = (t(:,cell_ii) - m)/std_t;
                        
                        m = mean(tt(start_frame:last_stim_frame,cell_ii),'omitnan');
                        std_t = std(tt(start_frame:last_stim_frame,cell_ii),'omitnan');
                        tt(:,cell_ii) = (tt(:,cell_ii) - m)/std_t;
                    end
                    %                 tt = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f(:,cell_IDs_comb{cell_i});
                    cell_count = cell_count+1;
                    All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common(:,cell_count) = mean(t,2,'omitnan');
                    All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common(:,cell_count) = mean(tt,2,'omitnan');
                    cells_done = [cells_done,cell_IDs_comb{cell_i}];
                end
            end
%             sum(cells_done<22)
%             size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common,2)
        end
    end
end
%%
save(fullfile(fid,'All_data.mat'),'All_data2', '-v7.3');
%% make a waterfall
x = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.Time_x;
y = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.cell_f_filt_common(:,1:50);
x = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.Time_x;
y = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.cell_f_filt_common(:,[3:3:152]);
waterfallplot(x,y,1,21);

% All_data = All_data2;
% All_data2 = table2cell(All_data2);
% save(fullfile(fid,'All_data_cell.mat'),'All_data2', '-v7.3');
% All_data2 = All_data;
%% Check stim allignment
for var_i = 1:length(varnames)
    cell_f_mat = nan(dims(1),dims(2),max_cell_no,max_frame_no);
    for animal_i = 1%:4%height(All_data2)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                figure; plot(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common)
                for stim_i = 1:length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use)
                    stimx = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use(stim_i);
                    xline(stimx);
                end
            end
        end
    end
end

%% Look at stims
frames_around_stim = 100;
pre_rmat = [];
post_rmat = [];
vars_mat = {};
for var_i = 1:length(varnames)
    for animal_i = 5:height(All_data2)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                sorted_stimframes = round(sort(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use)/6);
                t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common;
%                 t = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common;
                Pre_IP_stim_frames = sorted_stimframes(find(sorted_stimframes<mean(sorted_stimframes)));
                for stim_i = 1:length(Pre_IP_stim_frames)
                    if Pre_IP_stim_frames(stim_i)-frames_around_stim>0 && Pre_IP_stim_frames(stim_i)+frames_around_stim<length(t)
                        pre_rmat = [pre_rmat,t(Pre_IP_stim_frames(stim_i)-frames_around_stim:Pre_IP_stim_frames(stim_i)+frames_around_stim,:)];
                    end
                end
                
                Post_IP_stim_frames = sorted_stimframes(find(sorted_stimframes>mean(sorted_stimframes)));
                for stim_i = 1:length(Post_IP_stim_frames)
                    if Post_IP_stim_frames(stim_i)-frames_around_stim>0 && Post_IP_stim_frames(stim_i)+frames_around_stim<length(t)
                        post_rmat = [post_rmat,t(Post_IP_stim_frames(stim_i)-frames_around_stim:Post_IP_stim_frames(stim_i)+frames_around_stim,:)];
                    end
                end
            end
        end
    end
    vars_mat{var_i,1} = pre_rmat;
    vars_mat{var_i,2} = post_rmat;
end
% %%
figure;
count = 1;
x = linspace(1/5.2,201/5.2,201)-100/5.2;
for var_i = 1:size(vars_mat,1)
    for cond_i = 1:size(vars_mat,2)
        m = mean(vars_mat{var_i,cond_i},2,'omitnan');
        e = nansem_nansuite(vars_mat{var_i,cond_i},2);
        subplot(2,2,count); count = count +1;
        shadedErrorBar(x,m,e,'lineprops', {'color', params.colour_list{count}},'patchSaturation',0.1); hold on;
    end
end
%%
%             cell_ID_list = [];
%             plane_list = [];
%             ROI_list = [];
%             %             cell_ID_list = All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_ID]
%             for plane_i = 1:5
%                 cell_ID_list = [cell_ID_list,All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID'];
%                 plane_list = [plane_list,ones(1,length(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID))*plane_i];
%                 ROI_list = [ROI_list,linspace(1,length(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID),length(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID))];
%             end
%             cell_ID_mat =[cell_ID_mat;cell_ID_list];
%             plane_mat = [plane_mat;plane_list];
%             ROI_mat = [ROI_mat;plane_list];
%             var_exp_used{end+1,1} = varnames{var_i};
%             var_exp_used{end,2} = exp_list{exp_i};
            
%         end
%     end
    common_IDs = [];
%     count = 1;
    for cell_i = 1:size(cell_ID_mat,2)
        if length(unique(cell_ID_mat(:,cell_i))) == 1
            common_IDs(cell_i) = cell_ID_mat(1,cell_i);
        else
            cell_ID_mat(:,cell_i)
            % check fluorescence STD and keep the ROI with higher
            unique_tags = unique(cell_ID_mat(:,cell_i))   
            c = arrayfun(@(x)length(find(cell_ID_mat(:,cell_i) == x)), unique(cell_ID_mat(:,cell_i)), 'Uniform', false);
            cell2mat(c);
            
            for tag_i = 1:length(unique_tags)
                selected_exps = find(cell_ID_mat(:,cell_i) == unique_tags(tag_i));
                
                cell_STD
                for exp_i = 1:size(var_exp_used,1)
                    All_data.(var_exp_used{exp_i,1}){animal_i}.(var_exp_used{exp_i,2})
                end
            end
%             plane_mat(:,cell_i)
            fdsf
        end
            
    end
% end

%%
max_frame_no = []; max_cell_no = [];
for animal_i = 1:height(All_data2)
    for var_i = 1:length(varnames)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                    max_frame_no(animal_i,exp_i) = size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common,1);
                    max_cell_no(animal_i,exp_i) = size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common,2);
            end
        end
    end
end
dims = size(max_frame_no);
max_frame_no = max(max(max(max_frame_no))); 
max_cell_no = max(max(max(max_cell_no)));


%% Cells
% cell_f_mat = nan(dims(1),dims(2),max_cell_no,max_frame_no);
res_struct = [];
for var_i = 1:length(varnames)
    cell_f_mat = nan(dims(1),dims(2),max_cell_no,max_frame_no);
    for animal_i = 4%:4%height(All_data2)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
%                 for plane_i = 1:5
                    temp =  All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common;
                    % preprocess_cells_here
                    
                    sorted_stimframes = round(sort(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).stim_frames_to_use)/6);
                    last_stim_frame = sorted_stimframes(length(sorted_stimframes)/2-1);
                    start_frame = 5*60* round(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).fps / 6);
                    post_IP_frame = sorted_stimframes(length(sorted_stimframes)/2+1);
%                     post_IP_frame = 1;
                    for cell_i = 1:size(temp,2)
                        m = mean(temp(start_frame:last_stim_frame,cell_i),'omitnan');
                        std_t = std(temp(start_frame:last_stim_frame,cell_i),'omitnan');
                        temp(:,cell_i) = (temp(:,cell_i) - m)/std_t;
                    end
                    temp = temp(post_IP_frame:end,:);
                    cell_f_mat(animal_i,exp_i,1:size(temp,2),1:size(temp,1)) = temp';
                    x = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).Time_x;
                    xkeep(1:length(x)) = x;
%                 end
            end
        end
    end
    res_struct.(varnames{var_i}) = cell_f_mat;
end
% %%
figure;
for var_i = 1:length(varnames)
%    permuted_data = data_to_permutedatafun(res_struct.(varnames{var_i}));
temp = squeeze(mean(res_struct.(varnames{var_i}),2,'omitnan'));

permuted_data = data_to_permutedatafun(temp);
   m = mean(permuted_data,1,'omitnan');
   e = nansem_nansuite(permuted_data,1);
%    x = linspace(1,length(m),length(m));
   shadedErrorBar(xkeep,m,e,'lineprops', {'color', params.colour_list{var_i}},'patchSaturation',0.1); hold on;
end

%% Running

max_run_frames = []; 
for animal_i = 1:height(All_data)
    for var_i = 1:length(varnames)
        if ~isempty(All_data.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                    max_run_frames(animal_i,exp_i,plane_i) = size(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).abs_run_vec,1);
            end
        end
    end
end

dims = size(max_run_frames);
max_frame_no = max(max(max(max_run_frames))); 

run_mat = nan(dims(1),dims(2),dims(3),max_frame_no);
res_struct = [];
for var_i = 1:length(varnames)
    run_mat = nan(dims(1),dims(2),dims(3),max_frame_no);
    for animal_i = 1:height(All_data)
        if ~isempty(All_data.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                    temp =  All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).abs_run_vec;
                    
                    % preprocess_cells_here
                    x = All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).running_X;
%                     idx1 = find(x< All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).IP_timing_minutes(1));
                    idx2 = find(x> All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).IP_timing_minutes(2));
                    
%                     last_stim_frame = sorted_stimframes(length(sorted_stimframes)/2-1);
%                     start_frame = 1*60* round(All_data.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).fps / 6);
%                     post_IP_frame = sorted_stimframes(length(sorted_stimframes)/2+1);
%                     for cell_i = 1:size(temp,2)
%                         m = mean(temp(start_frame:last_stim_frame,cell_i),'omitnan');
%                         std_t = std(temp(start_frame:last_stim_frame,cell_i),'omitnan');
%                         temp(:,cell_i) = (temp(:,cell_i) - m)/std_t;
%                     end
                    temp = temp(idx2);
                    run_mat(animal_i,exp_i,1:length(temp)) = temp;
                    
            end
        end
    end
    res_struct.(varnames{var_i}) = run_mat;
end

figure;
for var_i = 1:length(varnames)
   permuted_data = data_to_permutedatafun(res_struct.(varnames{var_i}));
   m = mean(permuted_data,1,'omitnan');
   e = nansem_nansuite(permuted_data,1);
   x = linspace(1,length(m),length(m));
   shadedErrorBar(x,m,e,'lineprops', {'color', params.colour_list{var_i}},'patchSaturation',0.1); hold on;
end



%%
function [permuted_data] = data_to_permutedatafun(data)
%%
        % concatenates dimensions
        permuted_data = nan(size(data,1)*size(data,2),size(data,3));
        count = 0;
        for dim_1 = 1:size(data,1)
            for dim_2 = 1:size(data,2)
%                 for dim_3 = 1:size(data,3)
%                     for dim_4 = 1:size(data,4)
                        count = count+1;
                        permuted_data(count,:) = (data(dim_1,dim_2,:));
%                     end
%                 end
            end
        end
end

function subFolderNames = getsubfolder_names_fun(fid)
    % Get a list of all files and folders in this folder.
    files = dir(fid);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags); % A structure with extra info.
    % Get only the folder names into a cell array.
    subFolderNames = {subFolders(3:end).name};
end
function waterfallplot(x,y,filt_order,filt_wdw)
    figure;
    if size(y,1)<size(y,2)
        y = y';
    end
    for i = 1:size(y,2)
        plot(x,zscore(sgolayfilt(y(:,i),filt_order,filt_wdw))-i*3); hold on;
    end
    
end
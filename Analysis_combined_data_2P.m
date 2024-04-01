% Get all 2P data. Plot stims and runs
%%%%%% Analysis workflow:
% cell_f_extraction_2P_doublecheckplanes_24012022a.m
% ExportPreprocesseddata_2PGlusensstim08022022.m
% Combine_2Pdata0802.m
% Analysis_combined_data_2P.m
% Locate2Pcells.m
% Chech_running_per_gluClasses.m

% Combine_2Pdata0802.m
clear all; clc; close all

dir_t = 'P:\Alexander\Summer 2019\To_do\Sensory_stim_exp';
fid_t = 'All_data.mat';
animal_list = {'077','082'};
animal_list2 = {'126b','152f','G35','G36','m143a'};
load(fullfile(dir_t,fid_t));
colour_list = {[0 0 0 ],[ 1 0.07 0.07],[0.05 0.25 1], [1,0.72,0],[ 0.4940 0.1840 0.5560],[ 0.4660 0.6740 0.1880],[ 0.3010 0.7450 0.9330],[ 0.6350 0.0780 0.1840],[0,0,0],[0.5,0,0],[0,0.5,0],[0,0,0.5]};
    
%%
%% make a waterfall
x = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.Time_x;
y = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.cell_f_filt_common(:,1:50);
x = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.Time_x;
y = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.cell_f_filt_common(:,[3:3:152]);
waterfallplot(x,y,1,21);

% All_data2 = All_data2;
% All_data2 = table2cell(All_data2);
% save(fullfile(fid,'All_data_cell.mat'),'All_data2', '-v7.3');
% All_data2 = All_data2;
%% Check stim allignment
varnames = fieldnames(All_data2);
varnames = varnames(2:3);
max_run_frames = []; 
for animal_i = 1:height(All_data2)
    for var_i = 1:length(varnames)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                for plane_i = 1:5
                    max_run_frames(animal_i,exp_i,plane_i) = size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).abs_run_vec,1);
                    max_cell_no(animal_i,exp_i) = size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common,2);
                end
            end
        end
    end
end
max_cell_no = [];
dims = size(max_run_frames);
max_frame_no = max(max(max(max_run_frames))); 
max_cell_no = max(max(max(max_cell_no)));
%% plot mean HON signal around stim per experiment
frames_around_stim = [50,50]; % frames
runSamples_around_stim = [100,100]; % samples
idxs = [1:60;61:120];
clear Stim_table
Stim_table.ID = All_data2.ID;

for var_i = 1:length(varnames)
    cell_f_mat = nan(dims(1),dims(2),max_cell_no,max_frame_no);
    for animal_i = 1:height(All_data2) 
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                figure; 
                subplot(3,1,1);
                Time_x = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).Time_x * 60;
                
                plot(Time_x,mean(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_common,2)); hold on;
                stim_mat_sig = nan(frames_around_stim(1)+1+frames_around_stim(2),length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use),size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common,2));
                stim_mat_run = nan(runSamples_around_stim(1)+1+runSamples_around_stim(2),length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use));
             
                run_abs = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).abs_run_vec;
                runtime = linspace(0,Time_x(end),length(run_abs));
                ylabel('HON');
                yyaxis right;
                plot(runtime,run_abs); xlim([0, runtime(end)]);
                
                for stim_i = 1:length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use)                    
                    stimx = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_frames_to_use(stim_i);
                    stimx = round(stimx/6);
                    frames_to_take = [stimx-frames_around_stim(1):stimx+frames_around_stim(2)]; 
                    if max(frames_to_take)<= length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common)
                        stim_mat_sig(:,stim_i,:) = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_f_filt_common(frames_to_take,:);
                        
                        
                        stim_t = min(find(runtime>Time_x(stimx)));
                        
                        
                        frames_to_take = [stim_t-runSamples_around_stim(1):stim_t+runSamples_around_stim(2)];
                        stim_mat_run(:,stim_i) = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).abs_run_vec(frames_to_take);
                        
                        xline(Time_x(stimx));
                        
                    end
                end
                
                ylabel('Run');
                xlabel('Time (s)');
                
                
                subplot(3,1,2);
                fpst = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).fps/6;
                x = linspace(-frames_around_stim(1)/fpst,frames_around_stim(2)/fpst,size(stim_mat_sig,1));
                x_sig = x;
                mean_t = mean(mean(stim_mat_sig(:,idxs(1,:),:),3,'omitnan'),2,'omitnan');
                sem_t = nansem_nansuite(mean(stim_mat_sig(:,idxs(1,:),:),3,'omitnan'),2);
                shadedErrorBar(x,mean_t,sem_t,'lineprops', {'color', colour_list{1}},'patchSaturation',0.1); hold on;
                mean_t = mean(mean(stim_mat_sig(:,idxs(2,:),:),3,'omitnan'),2,'omitnan');
                sem_t = nansem_nansuite(mean(stim_mat_sig(:,idxs(2,:),:),3,'omitnan'),2);
                shadedErrorBar(x,mean_t,sem_t,'lineprops', {'color', colour_list{2}},'patchSaturation',0.1); hold on;
                xline(0);
                
                ylabel('HONs');
                xlabel('Time (s)');
                subplot(3,1,3);
                x = linspace(0,runtime(runSamples_around_stim(1)+runSamples_around_stim(2)+1),runSamples_around_stim(1)+runSamples_around_stim(2)+1);
                x = x-(x(end)/2);
                x_run = x;
                mean_t = mean(stim_mat_run(:,idxs(1,:),:),2);
                sem_t = nansem_nansuite(stim_mat_run(:,idxs(1,:),:),2);
                shadedErrorBar(x,mean_t,sem_t,'lineprops', {'color', colour_list{1}},'patchSaturation',0.1); hold on;
                mean_t = mean(stim_mat_run(:,idxs(2,:),:),2);
                sem_t = nansem_nansuite(stim_mat_run(:,idxs(2,:)),2);
                shadedErrorBar(x,mean_t,sem_t,'lineprops', {'color', colour_list{2}},'patchSaturation',0.1); hold on;
                xline(0);
                
                ylabel('Run speed');
                xlabel('Time (s)');
                
                name = strjoin({All_data2.ID{animal_i},varnames{var_i},exp_list{exp_i}});
                suptitle(strrep(name,'_',' '));
                
                Stim_table.(varnames{var_i}){animal_i}.(exp_list{exp_i}).Parameters = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).Parameters;
                Stim_table.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_mat_sig = stim_mat_sig;
                Stim_table.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_mat_sig = stim_mat_run;
                Stim_table.(varnames{var_i}){animal_i}.(exp_list{exp_i}).x_sig = x_sig;
                Stim_table.(varnames{var_i}){animal_i}.(exp_list{exp_i}).stim_mat_sig = stim_mat_run;
                
                yourFolder = fullfile('P:\Paulius\Glucose paper\2P analysis','StimFigs');
                if ~exist(yourFolder, 'dir')
                    mkdir(yourFolder);
                end
                filename = fullfile(yourFolder,name);
                saveas(gcf,filename);
                saveas(gcf,filename, 'png');
                close all
            end
        end
    end
end
safasfd
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
%             %             cell_ID_list = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i}).cell_ID]
%             for plane_i = 1:5
%                 cell_ID_list = [cell_ID_list,All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID'];
%                 plane_list = [plane_list,ones(1,length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID))*plane_i];
%                 ROI_list = [ROI_list,linspace(1,length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID),length(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(plane_i).cell_ID))];
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
                    All_data2.(var_exp_used{exp_i,1}){animal_i}.(var_exp_used{exp_i,2})
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
for animal_i = 1:height(All_data2)
    for var_i = 1:length(varnames)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                    max_run_frames(animal_i,exp_i) = size(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).abs_run_vec,1);
            end
        end
    end
end

dims = size(max_run_frames);
max_frame_no = max(max(max(max_run_frames))); 
plotting_timing = [-20,25];
running_speed_limit = 100;
movmean_vec = [0,10*60*3] ; % run fps is 10

run_mat = nan(dims(1),dims(2),14400);
res_struct = [];
for var_i = 1:length(varnames)
    run_mat = nan(dims(1),dims(2),max_frame_no);
    runx_mat = run_mat;
    for animal_i = 1:height(All_data2)
        if ~isempty(All_data2.(varnames{var_i}){animal_i})
            exp_list = fieldnames(All_data2.(varnames{var_i}){animal_i});
            for exp_i = 1:length(exp_list)
                    temp =  All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).abs_run_vec;
                    
                    % preprocess_cells_here
                    x = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).running_X;
                    %                     idx1 = find(x< All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).IP_timing_minutes(1));
                    if any(strcmpi(animal_list2,All_data2.ID{animal_i}))
                        IP_timing = All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).IP_timing_minutes(2) - 3;
                    else
                        IP_timing =  All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).IP_timing_minutes(1);
                    end
                    idx2 = find(x> IP_timing + plotting_timing(1) & x< IP_timing + plotting_timing(2)  );
%                     
%                     last_stim_frame = sorted_stimframes(length(sorted_stimframes)/2-1);
%                     start_frame = 1*60* round(All_data2.(varnames{var_i}){animal_i}.(exp_list{exp_i})(1).fps / 6);
%                     post_IP_frame = sorted_stimframes(length(sorted_stimframes)/2+1);
%                     for cell_i = 1:size(temp,2)
%                         m = mean(temp(start_frame:last_stim_frame,cell_i),'omitnan');
%                         std_t = std(temp(start_frame:last_stim_frame,cell_i),'omitnan');
%                         temp(:,cell_i) = (temp(:,cell_i) - m)/std_t;
%                     end
                    temp = temp(idx2);
                    run_mat(animal_i,exp_i,1:length(temp)) = temp;
                    runx_mat(animal_i,exp_i,1:length(temp)) = temp;
            end
        end
    end

    res_struct.(varnames{var_i}) = run_mat;
    res_structX.(varnames{var_i}) = runx_mat;
end
% %% 
close all
figure;

for animal_i  = 1:size(res_struct.(varnames{var_i}),1)
    for var_i = 1:length(varnames)
        subplot(size(res_struct.(varnames{var_i}),1),1,animal_i);
        temp = squeeze(res_struct.(varnames{var_i})(animal_i,:,:));
%         temp = squeeze(mean(res_struct.(varnames{var_i})(animal_i,:,:,:),2,'omitnan'));
        %    permuted_data = data_to_permutedatafun(temp);
        
        %    x_t = squeeze(mean(res_structX.(varnames{var_i}),2,'omitnan'));
        %    x = mean(data_to_permutedatafun(x_t),1,'omitnan');
        
        temp(:,find(isnan(mean(temp,1,'omitnan')))) = [];
        
        % remove weirdly high signal
        for exp_i = 1:size(temp,1)
            idx_to_nan = find(temp(exp_i,:)>running_speed_limit);
            if ~isempty(idx_to_nan)
                temp(exp_i,idx_to_nan) = nan;
            end
        end
        
        temp = movmean(temp,movmean_vec,2);
        m = mean(temp,1,'omitnan');
        if length(m>1)
            e = nansem_nansuite(temp,1);
            x = linspace( plotting_timing(1), plotting_timing(2),length(m));
            shadedErrorBar(x,m,e,'lineprops', {'color', colour_list{var_i}},'patchSaturation',0.1); hold on;
            xline(0);
        end
    end
end

figure;

for var_i = 1:length(varnames)
    
    temp = squeeze(res_struct.(varnames{var_i})(:,:,:));
    %         temp = squeeze(mean(res_struct.(varnames{var_i})(animal_i,:,:,:),2,'omitnan'));
    temp = data_to_permutedatafun(temp);
    
    %    x_t = squeeze(mean(res_structX.(varnames{var_i}),2,'omitnan'));
    %    x = mean(data_to_permutedatafun(x_t),1,'omitnan');
    
    temp(:,find(isnan(mean(temp,1,'omitnan')))) = [];
    
    % remove weirdly high signal
    for exp_i = 1:size(temp,1)
        idx_to_nan = find(temp(exp_i,:)>running_speed_limit);
        if ~isempty(idx_to_nan)
            temp(exp_i,idx_to_nan) = nan;
        end
    end
    
    temp = movmean(temp,movmean_vec,2);
    m = mean(temp,1,'omitnan');
    e = nansem_nansuite(temp,1);
    x = linspace( plotting_timing(1), plotting_timing(2),length(m));
    shadedErrorBar(x,m,e,'lineprops', {'color', colour_list{var_i}},'patchSaturation',0.1); hold on;
    xline(0);
end

%% make a waterfall
x = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.Time_x;
y = All_data2.Glu{2,1}.exp_20190827_13_50_09_082_stims_glu_IP_2708.cell_f_filt_common(:,1:50);
x = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.Time_x;
y = All_data2.Glu{1,1}.exp_20190828_16_52_52_077_stims_glu_IP_2808.cell_f_filt_common(:,[3:3:152]);
waterfallplot(x,y,1,21);


%%

function waterfallplot(x,y,filt_order,filt_wdw)
    figure;
    if size(y,1)<size(y,2)
        y = y';
    end
    for i = 1:size(y,2)
        plot(x,zscore(sgolayfilt(y(:,i),filt_order,filt_wdw))-i*3); hold on;
    end
    
end

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

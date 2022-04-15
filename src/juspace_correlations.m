
%% Calculate spatial correlation between a set of input and target maps 
% Spatial autocorrelation-corrected null-maps are created for *input_vols*

function [save_path] = juspace_correlations(input_vols, target_vols, atlas, n_perm, save_path, analysis)

% juspace options
if ~exist('analysis', 'var')
    analysis = 'spearman';
end
switch analysis
    case 'spearman'
        options = [4 1 NaN 4];
        disp('JuSpace: spearman correlations.')
    case 'mlr'
        options = [4 3 NaN 4];
        disp('JuSpace: multiple linear regression.')
    otherwise
        error('Define analysis type!')
end

% create temp folder
temp = fullfile(pwd, 'matlab_temp');
if ~exist(temp, 'dir')
    mkdir(temp);
end

% unzip function
function [unzipped_vols] = unzip_vols(vols)
    if isstring(vols) || ischar(vols)
        vols = cellstr(vols)
    end
    for i=1:length(vols) 
        [~, ~, ext] = fileparts(vols(i));
        if strcmp(ext, '.gz')
            disp(['Unzipping: ' vols{i}])
            unzipped_vols(i) = gunzip(vols(i), temp);
        else
            unzipped_vols(i) = vols(i);
        end
    end
end

% unzip volumes
disp('Unzipping input volume(s) (if necessary).')
input_vols = unzip_vols(input_vols);
disp('Unzipping target volume(s) (if necessary).')
target_vols = unzip_vols(target_vols);
disp('Unzipping atlas volume (if necessary).')
atlas = unzip_vols(atlas);

% JuSpace: Get data and correlate
[res,p_all,stats,data,D1,D2,target_data,Resh,T1] = compute_DomainGauges(input_vols',[],target_vols',atlas{1},options);

switch analysis
    case 'spearman'
        % JuSpace: Compute exact p value
        [p_exact,dist_rand] = compute_exact_spatial_pvalue(D1,target_data,atlas{1},res,n_perm,options,target_vols',T1);
        % results
        Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
        % make table
        T = cell2table(Resh(2:end,:), 'VariableNames', Resh(1,:));
    case 'mlr'
        % no permutation, make table
        T = cell2table(Resh(2:end,1:4), 'VariableNames', Resh(1,1:4));
end

% remove temp folder with unzipped volumes
rmdir(temp, 's')
disp('Deleted temp data.')

% save result table
writetable(T, save_path)

end

%% Download the data if necessary.
if ~exist('GSE44772.txt','file')
    % Only download data from net if file not already in directory.
    disp('Retrieving GSE data...')
    gsedata=getgeodata('GSE44772', 'ToFile','GSE44772.txt');
    disp('done!')
elseif ~exist('gsedata', 'var')
    % Only load from file if not already in workspace.
    disp('Loading GSE data...')
    gsedata=geoseriesread('GSE44772.txt');
    disp('done!')
else
    disp('Skipping download, reload of gene expression data')
end

if ~exist('genes.txt','file')
    % Only download data from net if file not already in directory.
    disp('Retrieving platform data...')
    genes=getgeodata(gsedata.Header.Series.platform_id,'ToFile','genes.txt');
    disp('done!')
elseif ~exist('genes', 'var')
    % Only load from file if not already in workspace.
    disp('Loading platform data...')
    genes=geosoftread('genes.txt');
    disp('done!')
else
    disp('Skipping download, reload of platform data')
end

%% filter data
expvalues=gsedata.Data;
gene=genes.Data;
[mask,Fdata] = genelowvalfilter(expvalues,'absval',log2(2));
[Fmask,Filtdata]=geneentropyfilter(Fdata,'Percentile',30);

%% Extract metadata

% Extract subject ids and tissue type.
reg_title = '(\d*)_(.*)';
subj_id = extract_meta(gsedata.Header.Samples.title, reg_title, 1);
tissue_type = extract_meta(gsedata.Header.Samples.title, reg_title, 2);

% Extract disease state.
reg_disease = 'disease status: (.*)';
disease_state = extract_meta(...
    gsedata.Header.Samples.characteristics_ch2(4,:), reg_disease, 1);

%% Find unique subject ids.
subj_id_unique = unique(subj_id);
[n_gene, ~] = size(Filtdata);

%% Stack data from 3 brain regions.
stacked_data = zeros(3*n_gene, numel(subj_id_unique));
for i = 1:numel(subj_id_unique)
    datarows = Filtdata(:, strcmp(subj_id, subj_id_unique(i)));
    for j = 1:3
        stacked_data(((j-1)*n_gene+1):j*n_gene) = ...
            datarows(:, j);
    end
end

%% heirarchical clustering
% use clustergram

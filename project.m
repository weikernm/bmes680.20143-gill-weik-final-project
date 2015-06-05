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
expvalues = gsedata.Data;
gene = genes.Data;
[mask,Fdata] = genelowvalfilter(expvalues,'absval',log2(2));
gene = gene(mask, :);
[Fmask,fildata] = geneentropyfilter(Fdata,'Percentile',30);
gene = gene(Fmask, :);
[h,p] = ttest(fildata');
i_sigp = p<0.001;
Filtdata = fildata(i_sigp,:);
gene = gene(i_sigp, :);
%i_nonan = sum(isnan(Filtdata), 1) == 0;
%Filtdata = Filtdata(i_nonan,:);
%gene = gene(i_nonan, :);
[~, n_genetypes] = size(gene);

%% Extract metadata

% Extract subject ids and tissue type.
reg_title = '(\d*)_(.*)';
subj_id = extract_meta(gsedata.Header.Samples.title, reg_title, 1);
tissue_type = extract_meta(gsedata.Header.Samples.title, reg_title, 2);

% Extract disease state.
reg_disease = 'disease status: (.*)';
disease_state = extract_meta(...
    gsedata.Header.Samples.characteristics_ch2(4,:), reg_disease, 1);
% Extract age.
reg_age = 'age: (\d*)';
age = str2double(extract_meta(...
    gsedata.Header.Samples.characteristics_ch2(1,:), ...
    reg_age, 1));
% Extract gender.
reg_gender = 'gender: ([MF])';
gender = extract_meta( gsedata.Header.Samples.characteristics_ch2(3,:), ...
    reg_gender, 1);

%% Find unique subject ids.
subj_id_unique = unique(subj_id);
[n_gene, ~] = size(Filtdata);

%% Stack data from 3 brain regions.
stacked_data = zeros(3*n_gene, numel(subj_id_unique));
stacked_genes = cell(3*n_gene, n_genetypes);
stacked_tissue = cell(3*n_gene, 1);
tissue_type_unique = unique(tissue_type);
for j = 1:3
     stacked_tissue(((j-1)*n_gene+1):j*n_gene) = tissue_type_unique(j);
     stacked_genes(((j-1)*n_gene+1):j*n_gene, :) = gene;
end
for i = 1:numel(subj_id_unique)
    datarows = Filtdata(:, strcmp(subj_id, subj_id_unique(i)));
    for j = 1:3
        stacked_data(((j-1)*n_gene+1):j*n_gene, i) = ...
            datarows(:, j);
    end
end

%% k-means clustering
opts = statset('Display','final');
[coeff, score, latent,~,explained,~]=pca(stacked_data');
[idx,ctrs] = kmeans(stacked_data',2,'Options',opts);
ctrsindx=ctrs*coeff;
figure
plot(score(idx==1,1),score(idx==1,2),'r.','MarkerSize',12)
hold on
plot(score(idx==2,1),score(idx==2,2),'b.','MarkerSize',12)
%hold on
%plot(score(idx==3,1),score(idx==3,2),'k.','MarkerSize',12)
hold on
plot(ctrsindx(:,1),ctrsindx(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrsindx(:,1),ctrsindx(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
hold off
title('K-means Clustering of Samples');
%% Cluster by disease and age
for i = 1:length(subj_id_unique)
   unique_id(i) = str2num(cell2mat(subj_id_unique(i)));
end
grps=disease_state(unique_id);
disease_state=strcmp(grps,'normal');
figure
plot(score(disease_state==1,1),score(disease_state==1,2),'r.','MarkerSize',12)
hold on
plot(score(disease_state==0,1),score(disease_state==0,2),'b.','MarkerSize',12)
hold on
plot(ctrsindx(:,1),ctrsindx(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrsindx(:,1),ctrsindx(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend('Normal','Alzeihmer''s Disease','Centroids',...
       'Location','NW')
   title('Cluster by Disease State and PCA score')
age_unique=age(unique_id);
[decsage,decageindx]=sort(age_unique);
diseaseage=disease_state(decageindx);
colorVec=winter(5);
group=repmat([1 2 3 4 5], 46, 1);
scoreage=score(decageindx,:);
figure
for i=1:230
    if diseaseage(1,i)==1
norm=plot(scoreage(i,1),scoreage(i,2),'*','Color',colorVec(group(i),:),'MarkerSize',12);
    else 
alz=plot(scoreage(i,1),scoreage(i,2),'o','Color',colorVec(group(i),:),'MarkerSize',12);
    end
hold on
end
 plot(ctrsindx(:,1),ctrsindx(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrsindx(:,1),ctrsindx(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)
legend([norm,alz],'Normal','Alzheimer''s disease','Centroids',...
       'Location','NW')
   title('Correlation between Age and Disease')
    hold off;
 hold off;
%% Seperate data into diseased and control groups
disease_control_tree = fitctree(stacked_data',grps','MaxCat',2);
resuberror = resubLoss(disease_control_tree)
view(disease_control_tree,'Mode','graph')

%% Correlate age and genes
% Assumes that missing data are normal.
Alz_score = zeros(size(grps));
Norm_score = zeros(size(grps));
Alz_score(strcmp(grps, 'Alzheimer''s disease')) = 1;
Norm_score(strcmp(grps, 'normal')) = 1 ;
Norm_age=age_unique(Norm_score==1);
Alz_age=age_unique(Alz_score==1);
[n_feature,~]=size(stacked_genes);
R_alz = zeros(n_feature, 1);
R_norm = zeros(n_feature, 1);
R_normage = zeros(n_feature, 1);
R_alzage = zeros(n_feature, 1);
P_alz = zeros(n_feature, 1);
P_norm = zeros(n_feature, 1);
P_normage = zeros(n_feature, 1);
P_alzage = zeros(n_feature, 1);
for i = 1:n_feature
   [R_alz_i, P_alz_i] = corrcoef(Alz_score', stacked_data(i, :));
   [R_norm_i, P_norm_i] = corrcoef(Norm_score', stacked_data(i, :));
       [R_age_i, P_age_i] = corrcoef(age_unique(~isnan(age_unique))', ...
      stacked_data(i, ~isnan(age_unique)));
   R_alz(i) = R_alz_i(2);
   R_norm(i) = R_norm_i(2);
   P_alz(i) = P_alz_i(2);
   P_norm(i) = P_norm_i(2);
   P_age(i) = P_age_i(2);
end
Age_sigp= P_age<0.001;
Alz_sigp = P_alz<0.001;
Norm_sigp = P_norm<0.001;
Agerelated=find(Age_sigp'==1 & Norm_sigp==1);
Alzrelated=find(Age_sigp'==1 & Alz_sigp==1);
genediff=setdiff(Alzrelated,Agerelated);
Age_Alz_genes=stacked_genes((genediff),4);


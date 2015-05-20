gsedata=getgeodata('GSE44772', 'ToFile','GSE44772.txt');
genes=getgeodata(gsedata.Header.Series.platform_id,'ToFile','genes.txt');

%% filter data
expvalues=gsdata.Data;
gene=genes.Data;
[mask, expvalues, genesfilt] = genelowvalfilter(expvalues,gene,...
                                              'absval',log2(4));
numel(genesfilt)

%% heirarchical clustering
% use clustergram

if ~exist('GSE44772.txt','file')
    gsedata=getgeodata('GSE44772', 'ToFile','GSE44772.txt');
else
    gsedata=geoseriesread('GSE44772.txt');
end
if ~exist('genes.txt','file')
    genes=getgeodata(gsedata.Header.Series.platform_id,'ToFile','genes.txt');
else
    genes=geosoftread('genes.txt');
end

%% filter data
expvalues=gsedata.Data;
gene=genes.Data;
[mask,Fdata] = genelowvalfilter(expvalues,'absval',log2(2));
[Fmask,Filtdata]=geneentropyfilter(Fdata,'Percentile',30);
%% heirarchical clustering
% use clustergram

% load gene expression density data
[data, tp_names, gene_names, struct_names]  = load_data();
% select the time period: E11.5, E13.5, E15.5, E18.5, P4, P14, P28
P14=squeeze(data(6,:,:));
P4=squeeze(data(5,:,:));
P28=squeeze(data(7,:,:));
diff_2814=P28-P14;
diff=P14-P4;
diff(:,12)=1:2002';

diff_cleaned=rmmissing(diff);
gene_cleaned=gene_names(diff_cleaned(:,12));
load("mri4.mat");
load("mri14.mat");

%% calculate combined mean and std for 11 ROIs
grpn = d(:,1); % number of voxels
ADmean = d(:,3); % image mean for each parameter
ADstd = d(:,4); % image stdev for each parameter
FAmean = d(:,5); % image mean for each parameter
FAstd = d(:,6); % image stdev for each parameter
MDmean = d(:,7); % image mean for each parameter
MDstd = d(:,8); % image stdev for each parameter
RDmean = d(:,9); % image mean for each parameter
RDstd = d(:,10); % image stdev for each parameter
NDImean = d(:,11); % image mean for each parameter
NDIstd = d(:,12); % image stdev for each parameter
ODImean = d(:,13); % image mean for each parameter
ODIstd = d(:,14); % image stdev for each parameter

[~, ADallmean, ADallstd] = overallmeanstd(grpn, ADmean, ADstd);
[~, FAallmean, FAallstd] = overallmeanstd(grpn, FAmean, FAstd);
[~, MDallmean, MDallstd] = overallmeanstd(grpn, MDmean, MDstd);
[~, NDIallmean, NDIallstd] = overallmeanstd(grpn, NDImean, NDIstd);
[~, ODIallmean, ODIallstd] = overallmeanstd(grpn, ODImean, ODIstd);
[~, RDallmean, RDallstd] = overallmeanstd(grpn, RDmean, RDstd);
result=[ADallmean,FAallmean,MDallmean,RDallmean,NDIallmean,ODIallmean,ADallstd,...
    FAallstd,MDallstd,RDallstd,NDIallstd,ODIallstd];
%% correlation analysis
% sort by one brain region
diff_sorted=sortrows(diff,1);
gene_sorted=gene_names(diff_sorted(:,12));

ad_co=diff_cleaned(diff_cleaned(:,1)>0&diff_cleaned(:,2)<0....
    &diff_cleaned(:,3)>0&diff_cleaned(:,4)>0&diff_cleaned(:,5)>0....
    &diff_cleaned(:,6)>0&diff_cleaned(:,7)>0&diff(:,8)>0....
    &diff_cleaned(:,9)>0&diff_cleaned(:,10)>0&diff_cleaned(:,11)>0,:);
ad_co_name=gene_names(ad_co(:,12));

% calculate correlation coeffieent and pvalues of gene expression vs. one 
% MRI metrics
clear r p Index gene_selected gene_selected_n p_selected d rank result gene_sorted;
for i=1:2002
    [R,P]=corrcoef([mri_14(:,4);mri_4(:,4)],[P14(i,:)';P4(i,:)']);
    r(i,1)=R(1,2);
    p(i,1)=P(1,2);
end

% select significant correlated genes and genes with averaged abs expression difference>=exp(1)
d=diff(:,1:11);
mean_difference=mean(abs(d),2);
Index=(p<=0.05 & mean_difference>=2.72); 
gene_selected=gene_names(Index);
gene_selected_n=diff(Index,12);
p_selected=(p(Index));
r_selected=r(Index);

% sort by r value
rank=[gene_selected_n mean_difference(Index) p_selected r_selected];
result=sortrows(rank,4,'descend');
gene_sorted=gene_names(result(:,1)); 

C=intersect(A,B);

% scatterplot MRI metric VS. gene expression of one gene with linear
% regression
for i=1:9
    figure(i)
    c = linspace(1,10,11);
    scatter(mri_14(:,4),P14(in(i),:),80,c,'filled')
    hold on
    scatter(mri_4(:,4),P4(in(i),:),80,c)
    % add a linear regression line
    mdl = fitlm([mri_14(:,4);mri_4(:,4)],[P14(in(i),:)';P4(in(i),:)']);
    plot(mdl)
    hold off
    % legend('RSP-P14','Tel-P14','HYP-P14','P3-P14','P2-P14','P1-P41','M-P14',....
        % 'PPH-P14','PH-P14','PMH-P14','MH-P14','RSP-P4','Tel-P4','HYP-P4',....
        % 'P3-P4','P2-P4','P1-P4','M-P4','PPH-P4','PH-P4','PMH-P4','MH-P4')
    title('RD VS. Expression of ',char(gene_names(in(i))))
    saveas(gcf,['RD vs ',char(gene_names(in(i))),'.jpg'])
end

% plot correlation between MRI metrics and gene expression
for i=1:10
    figure(i)
    X=[[mri_14(:,4);mri_4(:,4)],[P14(in(i),:)';P4(in(i),:)']];
    [R,Pvalue] = corrplot(X,VarNames={'RD' char(gene_names(in(i)))});
    %title('RD vs',char(gene_names(in(i))))
    saveas(gcf,['correlation RD vs ',char(gene_names(in(i))),'.jpg'])
end

% heatmap of gene expression difference in p4/p14
figure(2)
subplot(1,2,1)
imagesc(diff(:,1:11))
title('P14-P4')
xlabel('regions')
ylabel('genes')
colorbar
colorbarpwn(-20,15)
subplot(1,2,2)
imagesc(diff_2814)
colorbar
colorbarpwn(-10,10)
title('P28-P14')
xlabel('regions')


% select genes which are consisitently increasing/decreasing along time
r1_in=diff(diff(:,2)>0&diff_2814(:,2)>0,12);
r1_in_data=[P4(r1_in,1),P14(r1_in,1),P28(r1_in,1)];
r1_in_name=gene_names(r1_in);
r1_de=diff(diff(:,2)<0&diff_2814(:,2)<0,12);
r1_de_data=[P4(r1_de,1),P14(r1_de,1),P28(r1_de,1)];
r1_de_name=gene_names(r1_de);
r1_in_data_sorted=sortrows([r1_in_data,r1_in],[1,2,3]);
imagesc(r1_in_data_sorted(:,1:3))
colorbar

%% bargraph of MRI metrics changes (grouped by ROI)
load("mri4e.mat");
load("mri14e.mat");
X = categorical({'RSP','Tel','HYP','P3','P2','P1','M','PPH','PH','PMH','MH'});
X = reordercats(X,{'RSP','Tel','HYP','P3','P2','P1','M','PPH','PH','PMH','MH'});
% data = [p4(11,:)' p14(11,:)'];
% err = [p4e(11,:)' p14e(11,:)'];
data = [mri_4(:,6) mri_14(:,6)];
err = [mri_4e(:,6) mri_14e(:,6)];

b=bar(X,data);                

hold on

[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',data,err,'k','linestyle','none');

hold off
title('ODI')
legend('P4','P14')
saveas(gcf,'P4P14 MRI comparison-ODI.jpg')
%% bargraph of MRI metrics changes (grouped by P4/P14)
load("mri4e.mat");
load("mri14e.mat");
X = categorical({'P4','P14'});
X = reordercats(X,{'P4','P14'});
% data = [p4(11,:)' p14(11,:)'];
% err = [p4e(11,:)' p14e(11,:)'];
data = [mri_4(:,1) mri_14(:,1)]';
err = [mri_4e(:,1) mri_14e(:,1)]';

b=bar(X,data);                
colororder([jet(11)])
hold on

[ngroups,nbars] = size(data);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',data,err,'k','linestyle','none');
%ylim([0 0.6])

hold off
title('AD')
%legend('RSP','Tel','HYP','P3','P2','P1','M','PPH','PH','PMH','MH')
% saveas(gcf,'P4P14 MRI comparison-ODI.jpg')

%% partial least square regression (gene VS alltogether or each DTI separately)
X=rmmissing([P14';P4'],2);
y=[mri_14;mri_4]; % alltogether
% y=[mri_14(:,4);mri_4(:,4)]; % each DTI separately
ncomp=10;
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X,y,ncomp);
%Calculate the normalized PLS weights
W0 = stats.W ./ sqrt(sum(stats.W.^2,1)); 
%Calculate the variable importance in projection (VIP) scores for ncomp components.
p = size(XL,1); 
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
%Find variables with a VIP score greater than or equal to 1.
indVIP = find(vipScore >= 1);
gene_VIP=gene_cleaned(indVIP);

% sort by VIP score
vipScore_sorted=sortrows([vipScore,diff_cleaned(:,12)],1,'descend');
gene_sorted=gene_names(vipScore_sorted(:,2));

% %Plot the VIP scores.
% scatter(1:length(vipScore),vipScore,'x')
% hold on
% scatter(indVIP,vipScore(indVIP),'rx')
% plot([1 length(vipScore)],[1 1],'--k')
% hold off
% axis tight
% xlabel('Predictor Variables')
% ylabel('VIP Scores')

%plot correlation plot of PL1/2 VS. DTI
tit = ''; % figure title
gnames = {'P14','P4'}; % names of groups in data {dimension 1 and 2}
corrinfo = {}; % stats to display of correlation scatter plot
for i=1:2
    % XX=[[mri_14(:,6);mri_4(:,6)],XS(:,i)];
    % [R,Pvalue] = corrplot(XX,VarNames={'ODI',num2str(i)});
    data1=[mri_14(:,2),mri_4(:,2)]; % 1:AD,2:FA,3:MD,4:RD,5:NDI,6:ODI
    data2=[XS(1:11,i),XS(12:22,i)];
    label = {'','PL',num2str(i)}; % Names of data sets
    [cr, fig(i), statsStruct] = correlationPlot(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'colors',[0 0 1;1 0 0],'markerSize',10,'axesLimits','auto','symbols','..');
    % correlationPlot(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'colors',[0 0 1;1 0 0],'markerSize',10,'axesLimits','auto','symbols','..');
    set(gca,'XLim',[0.1 0.7],'YLim',[-0.4 0.4])
    set(gca,'XTick',(0.1:0.05:0.7),'YTick',(-0.4:0.2:0.4))
    axis off
    % exportgraphics(gcf,['correlation FA vs PL',num2str(i),'.png'],'BackgroundColor','none','ContentType','vector')
    saveas(gcf,['correlation ODI vs PL',num2str(i),'.jpg'])
end
%get correlation r-value and p-value of PL1-10 VS. DTI
for i=1:10
    [R,P]=corrcoef([mri_14(:,6);mri_4(:,6)],XS(:,i));
    r(1,i)=R(1,2);
    p(1,i)=P(1,2);
end

%% PLS component map
% import ADMBA segmentation
V=niftiread('/Users/hanxiny/Desktop/Mouse_Brain_development paper/gene-mri/mri data/P04/Histology_P04/P4_annotation.nii');
I=squeeze(V(247,end:-1:1,:)); % sagittal slice #247 for P4, #250 for P14
J=imrotate(I,90); % rotate into original degree of view
figure(1)
imagesc(J);
colormap([1 1 1; jet(256)])
caxis([15570 17700])
title('P14 segmentation label')

C=unique(J(:)); % find all segmentation labels note: add column2 of C as PLS1, column3 of C as PLS2, column4 of C as ROI# from PLS component map.xlsx

% draw PLS component map
PLS1=zeros(size(J,1),size(J,2));
PLS2=zeros(size(J,1),size(J,2));
ROI=zeros(size(J,1),size(J,2));
for i=1:size(J,1)
    for j=1:size(J,2)
        PLS1(i,j)=C(C(:,1)==J(i,j),2);
        PLS2(i,j)=C(C(:,1)==J(i,j),3);
        ROI(i,j)=C(C(:,1)==J(i,j),4);
    end
end
figure(2)
h1=imagesc(PLS1);
set(h1, 'AlphaData', 1-isnan(PLS1))
title('PLS1 P4')
colorbar
colormap("gray")
caxis([-0.4 0.3])
axis off

figure(3)
h2=imagesc(PLS2);
set(h2, 'AlphaData', 1-isnan(PLS2))
title('PLS2 P14')
colorbar
colormap("gray")
caxis([-0.4 0.4])
axis off

figure(4)
imagesc(ROI);
colormap([1 1 1; jet(11)])
title('P4 ROI')
axis off
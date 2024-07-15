%% load gene expression density data
[data, tp_names, gene_names, struct_names]  = load_data();
% select the time period: E11.5, E13.5, E15.5, E18.5, P4, P14, P28
P14=squeeze(data(6,:,:));
P4=squeeze(data(5,:,:));
P28=squeeze(data(7,:,:));
diff_2814=P28-P14;
load("mri4.mat");
load("mri14.mat");

%% partial least square regression on P4/P14 separately
X14=rmmissing(P14',2);
X4=rmmissing(P4',2);
y14=mri_14; % alltogether
y4=mri_4; % alltogether
ncomp=10;
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X14,y14,ncomp);
%Calculate the normalized PLS weights
W0 = stats.W ./ sqrt(sum(stats.W.^2,1)); 
%Calculate the variable importance in projection (VIP) scores for ncomp components.
p = size(XL,1); 
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
%Find variables with a VIP score greater than or equal to 1.
indVIP = find(vipScore >= 1);
P4(:,12)=1:2002';
diff_cleaned=rmmissing(P4);
gene_cleaned=gene_names(P4(:,12));
gene_VIP=gene_cleaned(indVIP);

% sort by VIP score
vipScore_sorted=sortrows([vipScore,diff_cleaned(:,12)],1,'descend');
gene_sorted=gene_names(vipScore_sorted(:,2));

%% PLS plots-PLS1 and PLS2 in separate figures
%plot correlation plot of PL1/2 VS. DTI
tit = ''; % figure title
gnames = {'P4'}; % names of groups in data {dimension 1 and 2}
corrinfo = {}; % stats to display of correlation scatter plot
para=1;% 1:AD,2:FA,3:MD,4:RD,5:NDI,6:ODI,7:AK,8:KFA,9:MK,10:RK
for i=1:1
    data1=X4(:,para); 
    data2=XS(:,i);
    label = {'','PL',num2str(i)}; % Names of data sets
    [cr, fig(i), statsStruct] = correlationPlot(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'colors',[0 0 1],'markerSize',10,'axesLimits','auto','symbols','..');
    axis on
    exportgraphics(gcf,['correlation FA vs PL',num2str(i),'.png'],'BackgroundColor','none','ContentType','vector')
    saveas(gcf,['P14 correlation AD vs PL',num2str(i),'-2.jpg'])
end
%get correlation r-value and p-value of PL1-10 VS. DTI
for i=1:10
    [R,P]=corrcoef(X4(:,para),XS(:,i));
    r(1,i)=R(1,2);
    p(1,i)=P(1,2);
end
%% PLS component map
% import ADMBA segmentation
V=niftiread('/Volumes/Xinyue_Drive/backup_from_laptop/Mouse_Brain_development paper/gene-mri/mri data/P14/Histology_P14/P14_annotation.nii');
I=squeeze(V(250,end:-1:1,:)); % sagittal slice #247 for P4, #250 for P14
% V=niftiread('/Volumes/Xinyue_Drive/backup_from_laptop/Mouse_Brain_development paper/gene-mri/mri data/P04/Histology_P04/P4_annotation.nii');
% I=squeeze(V(247,end:-1:1,:)); % sagittal slice #247 for P4, #250 for P14
J=imrotate(I,90); % rotate into original degree of view
figure(1)
imagesc(J);
colormap([1 1 1; jet(256)])
caxis([15570 17700])
title('P4 segmentation label')

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
title('PLS1 P14')
colorbar
colormap("jet")
% caxis([-0.5 0.4])
axis off

figure(3)
h2=imagesc(PLS2);
set(h2, 'AlphaData', 1-isnan(PLS2))
title('PLS2 P14')
colorbar
colormap("jet")
% caxis([-0.9 0.5])
axis off

figure(4)
imagesc(ROI);
colormap([1 1 1; jet(11)])
title('P4 ROI')
axis off

%% validity of PLS residuals
error=MSE(:,2:11);
Parameter=10;
Res=stats.Yresiduals(:,Parameter);
figure(1)
stem(Res,'LineWidth',2,'MarkerSize',12);
saveas(gcf,'P4 PLS residuals-add X data-stem-RK.png')
figure(2)
hh=normplot(Res); %normality
hh(1).MarkerSize=12;
hh(1).LineWidth=2;
hh(2).LineWidth=2;
hh(3).LineWidth=2;
[T,P_homo,df] = BPtest([XS(:,1),mri_14(:,Parameter)]);
[H, P_norm, W] = swtest(Res);
saveas(gcf,'P4 PLS residuals-add X data-normplot-RK.png')

%% variance and MSE as a function of PLS components
figure(3)
yyaxis left
plot(1:10,cumsum(100*PCTVAR(2,:)),'-o','LineWidth',2,'MarkerSize',12);
xlabel('number of PLS components')
ylabel('Percent Variance Explained in y')
yyaxis right
plot(1:10,error(2,:),'-o','LineWidth',2,'MarkerSize',12)
ylabel('MSE')

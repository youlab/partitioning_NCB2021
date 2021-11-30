
data6 = zeros(6,6);
%%
data0_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6-0P');
data6(:,1) = reshape(data0_6_raw(15:16,1:3),6,1);
%%
data0002_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6-0.0002P');
data6(:,2) = reshape(data0002_6_raw(15:16,1:3),6,1);
%%
data001_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6-0.001P');
data6(:,3) = reshape(data001_6_raw(15:16,1:3),6,1);
%%
data005_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6.0.005P');
data6(:,4) = reshape(data005_6_raw(15:16,1:3),6,1);
%%
data02_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6-0.02P');
data6(:,5) = reshape(data02_6_raw(15:16,1:3),6,1);
%%
data1_6_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','6-0.1P');
data6(:,6) = reshape(data1_6_raw(15:16,1:3),6,1);

%%
data24 = zeros(24,6);
%%
data0_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0P');
data24(:,1) = reshape(data0_24_raw(15:18,1:6),24,1);
%%
data0002_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.0002P');
data24(:,2) = reshape(data0002_24_raw(15:18,1:6),24,1);
%%
data001_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.001P');
data24(:,3) = reshape(data001_24_raw(15:18,1:6),24,1);
%%
data005_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.005P');
data24(:,4) = reshape(data005_24_raw(15:18,1:6),24,1);
%%
data02_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.02P');
data24(:,5) = reshape(data02_24_raw(15:18,1:6),24,1);
%%
data1_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.1P');
data24(:,6) = reshape(data1_24_raw(15:18,1:6),24,1);

%%
data24 = zeros(24,6);
%%
data0_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0P');
data24(:,1) = reshape(data0_24_raw(15:18,1:6),24,1);
%%
data0002_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.0002P');
data24(:,2) = reshape(data0002_24_raw(15:18,1:6),24,1);
%%
data001_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.001P');
data24(:,3) = reshape(data001_24_raw(15:18,1:6),24,1);
%%
data005_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.005P');
data24(:,4) = reshape(data005_24_raw(15:18,1:6),24,1);
%%
data02_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.02P');
data24(:,5) = reshape(data02_24_raw(15:18,1:6),24,1);
%%
data1_24_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','24-0.1P');
data24(:,6) = reshape(data1_24_raw(15:18,1:6),24,1);

%%
data96 = zeros(96,6);
%%
data0_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0P');
data96(:,1) = reshape(data0_96_raw(15:22,1:12),96,1);
%%
data0002_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0.0002P');
data96(:,2) = reshape(data0002_96_raw(15:22,1:12),96,1);
%%
data001_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0.001P');
data96(:,3) = reshape(data001_96_raw(15:22,1:12),96,1);
%%
data005_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0.005P');
data96(:,4) = reshape(data005_96_raw(15:22,1:12),96,1);
%%
data02_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0.02P');
data96(:,5) = reshape(data02_96_raw(15:22,1:12),96,1);
%%
data1_96_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','96-0.1P');
data96(:,6) = reshape(data1_96_raw(15:22,1:12),96,1);

%%
data384 = zeros(384,6);
%%
data0_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0Pert');
data384(:,1) = reshape(data0_384_raw(15:30,1:24),384,1);
%%
data0002_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0.0002P');
data384(:,2) = reshape(data0002_384_raw(15:30,1:24),384,1);
%%
data001_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0.001P');
data384(:,3) = reshape(data001_384_raw(15:30,1:24),384,1);
%%
data005_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0.005P');
data384(:,4) = reshape(data005_384_raw(15:30,1:24),384,1);
%%
data02_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0.02P');
data384(:,5) = reshape(data02_384_raw(15:30,1:24),384,1);
%%
data1_384_raw =xlsread('02-23-20-Glucose+Segragation Gradient.xlsx','384-0.1P');
data384(:,6) = reshape(data1_384_raw(15:30,1:24),384,1);

%%
data1536 = zeros(1536,6);
%%
data0_1536_raw =xlsread('02-23-20-1536-0.xls');
data1536(:,1) = data0_1536_raw(:,6);
%%
data0002_1536_raw =xlsread('02-23-20-1536-0002.xls');
data1536(:,2) = data0002_1536_raw(:,6);
%%
data001_1536_raw =xlsread('02-23-20-1536-001.xls');
data1536(:,3) = data001_1536_raw(:,6);
%%
data005_1536_raw =xlsread('02-23-20-1536-005.xls');
data1536(:,4) = data005_1536_raw(:,6);
%%
data02_1536_raw =xlsread('02-23-20-1536-02.xls');
data1536(:,5) = data02_1536_raw(:,6);
%%
data1_1536_raw =xlsread('02-23-20-1536-1.xls');
data1536(:,6) = data1_1536_raw(:,6);

%%
data_all = {data6,data24,data96,data384,data1536};

%%
ind = 1;
figure(1)
for j = 1:6
    for i = 1:5
        subplot(6,5,ind)
        tmp = data_all{i}(:,j);
        histogram(tmp(:),'FaceColor',[1,1,1]*0.7)
        ind = ind+1;
        if i ==1
            axis([0.05,0.3,0,3])
        elseif i ==2
            axis([0.05,0.4,0,15])
        elseif i ==3
            axis([0.03,0.5,0,50])  
        elseif i ==4
            axis([0.03,0.6,0,200])    
        elseif i ==5
            axis([0.4,1.2,0,500])               
        end
    end
end

subplot(6,5,26)
xlabel('OD')
subplot(6,5,21)
ylabel('Count (# of wells)')

%% from following up experiments with 2 repeats to replace 3 panels
data1536 = zeros(1536,6);

% data0NEW_1536_raw =xlsread('1536_2-27-20-CA.xls');
data0NEW_1536_raw =xlsread('1536_2-27-20-CA0_2.xls');
data1536_new(:,1) = data0NEW_1536_raw(:,6);

data005NEW_1536_raw =xlsread('1536_2-27-20-CA0.005.xls');
data1536_new(:,4) = data005NEW_1536_raw(:,6);

data24NEW = zeros(24,6);

data001NEW_24_raw =xlsread('2-27-20.xlsx','24_0.001');
data24NEW(:,5) = reshape(data001NEW_24_raw(9:12,1:6),24,1);

%%
subplot(6,5,5)
histogram(data1536_new(:,1),'FaceColor',[1,1,1]*0.7)
axis([0.4,1.2,0,800]) 

subplot(6,5,12)
histogram(data24NEW(:,5),'FaceColor',[1,1,1]*0.7)
axis([0.05,0.4,0,10])

subplot(6,5,20)
histogram(data1536_new(:,4),'FaceColor',[1,1,1]*0.7)
axis([0.4,1.2,0,800])  

%%
figure(2)
imagesc((reshape(data1536_new(:,1),48,32)))

%%
% data1536 = zeros(1536,6);
% 

figure(4)
subplot(1,5,1)
data0NEW_6_raw =xlsread('2-27-20.xlsx','6_0');
hist(reshape(data0NEW_6_raw(9:10,1:3),6,1));
axis([0.05,0.3,0,3])

subplot(1,5,2)
data0NEW_24_raw =xlsread('2-27-20.xlsx','24_0');
hist(reshape(data0NEW_24_raw(9:12,1:6),24,1));
axis([0.05,0.4,0,10])
data24NEW = zeros(24,6);

subplot(1,5,3)
data0NEW_96_raw =xlsread('2-27-20.xlsx','96_0');
hist(reshape(data0NEW_96_raw(15:22,1:12),96,1))
axis([0.03,0.5,0,40]) 
 
subplot(1,5,4)
data0NEW_384_raw =xlsread('2-27-20.xlsx','384_0');
hist(reshape(data0NEW_384_raw(15:30,1:24),384,1))
axis([0.03,0.6,0,250])

subplot(1,5,5)
hist(data0NEW_1536_raw(:,6))
axis([0.4,1.2,0,800])  

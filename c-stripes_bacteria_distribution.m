% Written by Taras Turiv (tturiv@kent.edu, taras.turiv@gmail.com)

% ************************************************************
% The code is used to perform the analysis of spatial distribution of
% bacteria moving through the patterned director field of periodic
% 1D modulated splay-bend deformation
% ***********************************************************

%c-stripes (distribution of bacteria in 2D)
%scale = 4096/135.97/10; %p/um --> 20x (ELWD)
%% Concentration
scale = 8192/262.02/10*0.97; %px/um --> 20x (DF, 2.1mm WD)
h = 3100; w = 2300;
minArea = 10;
minAR = 3;
text_size = 18;
nFiles = 1200;
% data export
data = cell(nFiles,1);
for i = 1:nFiles
    data{i,1} = dlmread(['Cell_3_5_for_analysis/unmask_bright_spots/unmask_bright_spots/res_',num2str(i),'.csv'],',',1,2);
end
%% Concentration 
% Uniform c-stripes (Fig.1)
% Bacteria concetration distribution map
dStep = 5*scale; %floor(h/10);
H = 0:dStep:h;
nDist = length(H);
D = zeros(nDist,nFiles);
for i = 1:nFiles
    %if s.bytes < 1000
    %    D(:,k) = nan;
    %    continue
    %end   
    area = data{i,1}(:,1);
    AR = data{i,1}(:,18);
    X = data{i,1}(:,7); Y = data{i,1}(:,8);
    X = X(area>minArea & AR>minAR); %Y<maxY & Y>minY & 
    Y = Y(area>minArea & AR>minAR);
    nPoints = length(X);
    for j = 1:nDist-1
        D(j,i) = D(j,i) + length(find(Y >= H(j) & Y < H(j+1)));
    end
end
dV = dStep*w/scale^2 * 20 * 1e-18;
Time = [0,nFiles-1]/2;
Radius = [max(H), H(1)]/scale+44;
figure('Position',[680 530 600 420],'PaperPositionMode','auto');
image(Time,Radius,D/dV*1e-14,'CDataMapping','scaled'); 
set(gca,'YDir','normal','YTick',160:160:1000,'YLim',[115 1000],'YMinorTick','on','XMinorTick','on');
c = colorbar; c.Label.String = '{\itc}({\ity}) (\times10^{14} m^{-3})';
set(gcf,'Position',[680 530 560 420],'PaperPositionMode','auto','Color','w');
c.Label.FontSize = 14;
c.Label.Position = [-0.5 10];
c.Position = [0.9 0.25 0.03 0.3];
c.FontSize = 14;
set(gca,'FontSize',text_size,'LineWidth',1,'XTick',[0:200:600],'XColor','w','YColor','w','LineWidth',1.5); grid off; %,'YTickLabel',1:6
box on; %axis square; 
xlim([0 max(max(Time))]); %ylim([0 max(max(Radius))]); %ylim([0 1000]); %axis tight;
set(gca,'Units','pixels','Position',[125 100 350 350]);
set(gca,'XTickLabel',sprintfc('\\color{black}%g',get(gca,'XTick')),'YTickLabel',sprintfc('\\color{black}%g',get(gca,'YTick')));
set(gcf,'Position',[600 500 620 520],'PaperPositionMode','auto','Color','w');
xlabel('{\itt} (s)','Color','k');
ylabel(['{\ity} (',char(181),'m)'],'Color','k');
set(c,'Position',[0.8 0.22 0.03 0.3])
ylim([160+80 960-80]);
pause(2);
%% Time averaged concentration distribution 2D
figure('Position',[600 500 620 540],'PaperPositionMode','auto','Color','w');%smooth
errorbar(((max(H)-H)/scale)'+44,(mean(D,2)/dV*1e-14),std(D,0,2)/dV*1e-14,'o','MarkerFaceColor','k','MarkerEdgeColor','w','Color',[0 0 0],'MarkerSize',8,'LineWidth',1);
data1 = [((max(H)-H)/scale)'+44,(mean(D,2)/dV*1e-14)];
view([90 -90]);
hold on; 
%theoretical prediction
k=6; N=960; L=160*k; y = linspace(0,L*(1-1/N),N); Dc = 10; V0 = 10; alph=L*V0/(k*pi*Dc); 
cavg = 2.7688e+13;
cp=exp(-alph*cos(k*y*pi/L)); 
cm=exp(alph*cos(k*y*pi/L));

C=cavg*N/sum(cp+cm); %Scaling
cp=cp*C; cm=cm*C;
plot(linspace(0,2*960,960*2),[(cp+cm),flipud(cp+cm)]*1e-14,'LineWidth',2.5,'Color','r'); 
data2 = [linspace(0,2*960,960*2)',([(cp+cm),flipud(cp+cm)]*1e-14)'];

cp60=cp60(:,512); cm60=cm60(:,512); cp600=cp600(:,512); cm600=cm600(:,512);

%plot(0.27*cosh(160*10/(pi*100)*cos(pi*(x-2)/(160))),x,'r','LineWidth',1.5);
%plot((cm+cp)*1e-14,x,'r','LineWidth',2);
%cTot = cm600 + cp600;
%cTot = [cTot; flipud(cTot); cTot; flipud(cTot)];
%plot(cTot,(1:length(cTot))/(512/160),'LineWidth',2,'Color',[1 0.75 0]);
%set(gca,'Position',[0.13 0.1473 0.7750 0.7777],'XTick',0:1:5);

plot(linspace(1,960*2+1,1024*2),[(cp600+cm600)*cavg/0.05;flipud(cp600+cm600)*cavg/0.05]*1e-14,'LineWidth',2,'Color','b','LineStyle','--'); %[1 0.75 0] %plot(flipud(cp+cm)*1e-14,linspace(960,960+960,1024),'r','LineWidth',1.5);
plot(linspace(1,960*2+1,1024*2),[(cp60+cm60)*cavg/0.05;flipud(cp60+cm60)*cavg/0.05]*1e-14,'LineWidth',2,'Color',[0 1 0],'LineStyle','--');
plot(y1,(c1+c2)/1e14,'-'); hold on;

data3 = [linspace(1,960*2+1,1024*2)',([(cp600+cm600)*cavg/0.05;flipud(cp600+cm600)*cavg/0.05]*1e-14)];
data4 = [linspace(1,960*2+1,1024*2)',[(cp60+cm60)*cavg/0.05;flipud(cp60+cm60)*cavg/0.05]*1e-14];
data5 = [y1',((c1+c2)/1e14)'];

set(gca,'YTick',0:1:5,'XTick',0:160:1000,'YMinorTick','on','XMinorTick','on','FontSize',text_size);
%axis square; 
axis tight;
ylabel('\langle{\itc}({\ity})\rangle_{\itt} (\times10^{14} m^{-3})','FontSize',text_size,'Interpreter','tex'); %xlim([0 50])
%ylabel(['{\ity} ({\itL})'],'FontSize',text_size); zlabel('Count'); %,char(181),'m)'
xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size); zlabel('Count');
l = legend({'Experiment','\tau\rightarrow\infty','\tau=600 s','\tau=60 s','\gamma>0'},'FontSize',14,'Box','on','Position',[0.5 0.21 0.2635 0.15]);
%legend({'Experiment',['Dc=',num2str(Dc),' ',char(181),'m^2/s']},'FontSize',14,'Box','on','Position',[0.5651 0.2170 0.2635 0.1806]);
xlim([80 1020]);
ylim([0 6]);
set(gca,'Units','pixels','Position',[125 100 350 350]);
pause(2);
myLines = get(gca,'Children');
set(myLines(2:4),'LineWidth',1.5);
set(l,'Position',[0.625 0.21 0.25 0.15]);
hold on;
xlim([160+80 960-80]);
%% Long term bacteria
figure('Position',[600 500 620 540],'PaperPositionMode','auto','Color','w');%smooth
load('d1.mat');
semilogy(d1(:,1),d1(:,2),'o',...
    d1(:,1),d1(:,4)/2,'o',...
    d1(:,1),d1(:,3),'o','MarkerFaceColor','auto');
set(gca,'FontSize',18,'XTick',0:1800:100000,'XTickLabel',(0:1800:100000)/3600);
legend({'Length','Number of bacteria','Width'},'FontSize',14,'Location','northeast');
xlabel('{\itt} (hrs)');
xlim([0 15000]);
ylabel(['{\itn}, {\itw} (',char(181),'m), {\itl} (',char(181),'m)']);
%% data into csv
fname = 'file_data_4c.csv';
fid = fopen(fname,'w'); 
fprintf(fid,'y(um),c(y)(*10^(14) m^(-3))\n'); 
fclose(fid);
dlmwrite(fname,data4c,'-append');
%% Trajectories
figure;
text_size = 18;
[a,~] = size(TrackID);
%figure('Position',[200 1000 1200 800]);
for i=1:a
    if length(TrackID(i,1).col)>10 %&& max(TrackID(i,1).row/scale)>300 && max(TrackID(i,1).row/scale)<450 && max(TrackID(i,1).col/scale)>600
        x = (w-TrackID(i,1).col)/scale;
        y = (h-TrackID(i,1).row)/scale;
        p1 = [x(end-1),y(end-1)];
        p2 = [x(end),y(end)];
        dp = (p2-p1)/sqrt(sum((p2-p1).^2));
        clr = rand(1, 3);
        plot(x,y,'Color',clr,'LineWidth',2);%,'Tag','trajectories');
        if length(TrackID(i,1).col) > 50
            %quiver(p1(1), p1(2), dp(1), dp(2),'Color',clr,'MaxHeadSize',1);
        end
        hold on;
    end
end
axis equal;
%axis equal;
%axis tight;
set(gca,'FontSize',text_size,'XLim',[0,700],'YLim',[0,1000],...
    'YTick',(0:160:1000)-45,'YTickLabel',(0:160:1000),'YMinorTick','on','XMinorTick','on',...
    'XTick',(0:200:800));
ylabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
xlabel(['{\itx} (',char(181),'m)'],'FontSize',text_size);
%%
set(gca,'Units','pixels','Position',[75 100 420 420]);
set(gcf,'Position',[600 500 600 600],'PaperPositionMode','auto','Color','w');

%% Fig.2a add director field
scale = 8192/262.02/10*0.97; %px/um --> 20x (DF, 2.1mm WD)
%im = imread('/media/turbik/DATA/PhD/DATA/2017_Bacteria_c,s,v-stripes/NewRM_cells/Cell_3_5_for_analysis/Cell_3_5_for_analysis_20fps0000_inverted.jpg');
im = imread('/media/turbik/DATA/PhD/DATA/2017_Bacteria_c,s,v-stripes/NewRM_cells/Cell_3_5_for_analysis/Cell_3_5_for_analysis_20fps0000_inverted-1.bmp');
image(im); colormap(gray(256)); axis equal; axis tight; hold on;
[h,w,~] = size(im);
[x,y] = meshgrid(1:50:w,1:50:h);
nx = cos(pi*(y+250)/(160*scale*1.01));
ny = -sin(pi*(y+250)/(160*scale*1.01));
quiver(x,y,nx,ny,0.25,'Color',[1 0.75 0.75],'ShowArrowHead','off','LineWidth',2); hold on;
quiver(x,y,-nx,-ny,0.25,'Color',[1 0.75 0.75],'ShowArrowHead','off','LineWidth',2);
%ylim([0 2500]); xlim([0 1855]);
axis off;
%% Small trajectories
%POSITION_X(ID == i), POSITION_Y
scale = 1;
[a,~] = size(ID);
%figure('Position',[200 1000 1200 800]);
im = imread('/media/turbik/DATA/PhD/DATA/2017_Bacteria_c,s,v-stripes/NewRM_cells/Cell_3_5_for_analysis/im0001-1.jpg');
image(im); colormap(gray(256)); hold on;
index = zeros(1,1); k = 1;
l = [12,23,24,43,1409,1534,1568];
for i=1:a
    x = POSITION_X(ID == i);
    if length(x) > 10 
        x = ((x))/scale+10;
        y = ((POSITION_Y(ID == i)))/scale;
        %&& max(TrackID(i,1).row/scale)>300 && max(TrackID(i,1).row/scale)<450 && max(TrackID(i,1).col/scale)>600
        %x = (w-TrackID(i,1).col)/scale;
        %y = (h-TrackID(i,1).row)/scale;
        %x = (ID(i,1).col)/scale;
        %y = (ID(i,1).row)/scale;
        %p1 = [x(end-1),y(end-1)];
        %p2 = [x(end),y(end)];
        %dp = (p2-p1)/sqrt(sum((p2-p1).^2));
        
        %c = rand(1)*0.75; clr = [c,c,1];
        clr = rand(1,3);
        if i == 1409
            plot(smooth(x(1:8:end-71)),smooth(y(1:8:end-71)),'Color',clr,'LineWidth',0.5);%,'Tag','trajectories');
        elseif i == 23
            plot(smooth(x(1:8:end-220)),smooth(y(1:8:end-220)),'Color',clr,'LineWidth',0.5);%,'Tag','trajectories');
        else
            plot(smooth(x(1:8:end)),smooth(y(1:8:end)),'Color',clr,'LineWidth',0.5);%,'Tag','trajectories');
        end
%         if length(TrackID(i,1).col) > 50
%             %quiver(p1(1), p1(2), dp(1), dp(2),'Color',clr,'MaxHeadSize',1);
%         end
        
        hold on;
        index(k,1) = i;
        k = k + 1;
    end
end
legend(gca,{num2str(index)});
axis equal; axis tight;
%ylim([0 632]/scale); xlim([0 516]/scale);
%% Instant velocity on Fig.2B
scale = 8192/262.02/10*0.97;
filename = '/media/turbik/DATA/PhD/MATLAB/Bacteria/c-stripes/traject_1-3_frames1.csv';
startRow = 2;
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec,'Delimiter',',','HeaderLines',startRow-1);fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells
d = cell2mat(raw);
%clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

X = d(:,5); Y = d(:,6);
TRACK_ID = d(:,3);
maxTr = max(TRACK_ID);
minTr = min(TRACK_ID);
for i = minTr:maxTr
    xTrack = X(TRACK_ID == i);
    yTrack = Y(TRACK_ID == i);
    
    x0(i+1) = xTrack(1);
    y0(i+1) = yTrack(1);
    
%    plot(xTrack,yTrack); hold on;
    vx(i+1) = nanmean(diff(X(TRACK_ID == i)));
    vy(i+1) = nanmean(diff(Y(TRACK_ID == i)));   
end
%axis equal;

%quiver(x0,y0,vx,vy)
figure;
im = imread('Cell_3_5_for_analysis_20fps0000-2.png');
image(im); colormap(gray(256)); axis tight; axis equal; hold on;
[h,w,~] = size(im);
axis off; 
cond = abs(vx)>0;

cond = ...
(round(x0,1) == 159.6 & round(y0,1) == 131.6) |...
(round(x0,1) == 151.4 & round(y0,1) == 130.7) |...
(round(x0,1) == 210.7 & round(y0,1) == 111.3) |...
(round(x0,1) == 458.5 & round(y0,1) == 109.5) |...
(round(x0,1) == 626.4 & round(y0,1) == 96.0) |...
(round(x0,1) == 30.4  & round(y0,1) == 604.1) |...
(round(x0,1) == 160.5  & round(y0,1) == 605) |...
(round(x0,1) == 278.2  & round(y0,1) == 604) |...
(round(x0,1) == 410.6  & round(y0,1) == 603.1) |...
(round(x0,1) == 324.8  & round(y0,1) == 1074) |...
(round(x0) == 162  & round(y0) == 1091) |...
(round(x0) == 268  & round(y0) == 1084) |... %?
(round(x0) == 278  & round(y0) == 604) |...
(round(x0) == 51  & round(y0) == 1078) |...
(round(x0) == 556  & round(y0) == 108) |...
(round(x0) == 673  & round(y0) == 600) |...
(round(x0) == 585  & round(y0) == 1078);
%(round(x0) == 474  & round(y0) == 1086) |...
%(round(x0) == 544  & round(y0) == 594) |...

%annotation('arrow',[1-585/w .2],[1-1078/h .2]);
%(round(x0,1) == 673.1  & round(y0,1) == 599.9) |...
%(round(x0,1) == 509.6  & round(y0,1) == 597.6) |...

%cond = abs(vx(1:end-1))>5 & ((abs(diff(x0))>1 & y0(1:end-1)<300) | (abs(diff(x0))>100 & y0(1:end-1)>400 & y0(1:end-1)<800) | (abs(diff(x0))>100 & y0(1:end-1)>900)) & round(x0(1:end-1))~=543 & round(y0(1:end-1))~=612 & round(x0(1:end-1))~=411 & round(y0(1:end-1))~=603 & x0(1:end-1)<750 & round(x0(1:end-1))~=383 & round(y0(1:end-1))~=596; 
%quiver(x0(cond),y0(cond),vx(cond)/scale/0.5,vy(cond)/scale/0.5,0.5,'r','LineWidth',5,'AutoScale','on','MaxHeadSize',200)

[Ma,Imax] = max(vx);
ylim([500 1200]);
%quiver(x0(Imax),y0(Imax)-20,vx(Imax)/scale/0.5,vy(Imax)/scale/0.5,10,'g','LineWidth',1.5,'AutoScale','off')

%set(gca,'Units','pixels','Position',[0 0 972/1.5 1276/1.5]);
set(gca,'Units','normal','Position',[0 0 1 1]);
set(gcf,'Position',[0 0 972/1.5 1276/1.5],'PaperPositionMode','auto','Color','w');
%% Velocity
%(round(x0,1) == 151.4 & round(y0,1) == 130.7) |...
xx = [151.4,324.8,210.7,458.5,626.4];
yy = [130.7,1074,111.3,109.5,96.03];
for i = 1:length(xx)
    vxx(i) = vx(round(x0) == round(xx(i)) & round(y0) == round(yy(i)));
    vyy(i) = vy(round(x0) == round(xx(i)) & round(y0) == round(yy(i)));
    annotation('arrow',(xx(i)/w+0.04)*[1 1]+[0 vxx(i)/50],(1-yy(i)/h)*[1 1]+[0 -vyy(i)/50],'Color','b');
end
%% Velocity from PIV!!!
scale = 8192/262.02/10*0.97;
X = x; Y= y;
x = X{1,1};
y = Y{1,1};
vx = u_original{1,1};
vy = v_original{1,1};

[h,w] = size(x);

vxavgy = zeros(h,1);
vyavgy = zeros(h,1);

for i = 1:h
    xnavg = 0;
    ynavg = 0;
    for j = 1:w
        if ~isnan(vx(i,j))
            vxavgy(i) = vxavgy(i) + vx(i,j);
            xnavg = xnavg + 1;           
        end
        if ~isnan(vy(i,j))
            vyavgy(i) = vyavgy(i) + vy(i,j);
            ynavg = ynavg + 1;          
        end    
    end
    vxavgy(i) = vxavgy(i)/xnavg;
    vyavgy(i) = vyavgy(i)/ynavg;
end
% v average
[h,w] = size(X{1,1});
vx_avg = zeros(h,w);
vy_avg = vx_avg;

navg = vx_avg;

for k = 1:1199
    for i = 1:h
        for j = 1:w
            if ~isnan(u_original{k,1}(i,j)) && ~isnan(v_original{k,1}(i,j))
                vx_avg(i,j) = vx_avg(i,j) + u_original{k,1}(i,j);
                vy_avg(i,j) = vy_avg(i,j) + v_original{k,1}(i,j);
                navg(i,j) = navg(i,j) + 1;
            end    
        end
    end
end
vx_avg = vx_avg./navg;
vy_avg = vy_avg./navg;
% Limit the number of vectors with smaller statistics
vx_temp = vx_avg;
vy_temp = vy_avg;
statN = 15;
vx_temp(navg < statN) = nan;
vy_temp(navg < statN) = nan;
% 2D y axis distribution
vx_y = zeros(h,1);
vy_y = zeros(h,1);
nvy = zeros(h,1);
for i=1:h
    for j=100:200
        if ~isnan(vx_temp(i,j))
            vx_y(i,1) = vx_y(i,1) + vx_temp(i,j);
            vy_y(i,1) = vy_y(i,1) + vy_temp(i,j);
            nvy(i,1) = nvy(i,1) + 1;
        end
    end
end
yax = (h-1:-1:0)*10/scale; 
vx_y = vx_y./nvy*1.1906e-05/1.9725e-06;
vy_y = vy_y./nvy*1.1906e-05/1.9725e-06;
astart = 0;
%% v_x_t velocity plot
text_size = 18;
figure('Position',[680 530 600 400],'PaperPositionMode','auto');
% hold on;
plot(fliplr(yax(1:end-astart))+57.5,(vx_y(1:end-astart)/1e-6),'.-','LineWidth',1.5,'Color','b','MarkerSize',20);
data2g = [fliplr(yax(1:end-astart))'+57.5,(vx_y(1:end-astart)/1e-6),(vy_y(1:end-astart)/1e-6)];
hold on; 
plot(fliplr(yax(1:end-astart))+57.5,(vy_y(1:end-astart)/1e-6),'.-','LineWidth',1.5,'Color','r','MarkerSize',20);
%data2h = [fliplr(yax(1:end-astart))'+57.5,(vy_y(1:end-astart)/1e-6)];

%axis square;
axis tight; %grid on;
set(gca,'FontSize',text_size,'YLim',[-5 12],'XTick',0:160:1000,...
    'Unit','pixel','Position',[125 100 420 250],'XMinorTick','on','YMinorTick','on'); 
hold on; 

%plot(x1*160+160,y1*20,'-');

xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
ylabel(['\langle{\itv}_{{\itx},{\ity}}({\ity})\rangle_{\itt} (',char(181),'m/s)'],'FontSize',text_size);
l = legend({'{\itv}_{\itx}','{\itv}_{\ity}'},'Box','off','FontSize',14);%,'Position',[0.75    0.80    0.1262    0.2125]);
%set(gcf,'Position',[600 500 600 400],'PaperPositionMode','auto','Color','w');

%set(gca,'Units','pixels','Position',[125 100 420 420]);
set(gcf,'Position',[600 500 620 400],'PaperPositionMode','auto','Color','w');

%xlim([0 960+160]);

xlim([160+80 960-80]);
ylim([-5 15]);

plot(0:1200,zeros(1201,1),'k-');

pause(2); 
%% data write into csv
fid = fopen('file_data_2g.csv','w'); 
fprintf(fid,'y(um),v_x(y)(um/s),v_y(y)(um/s)\n');
fclose(fid);
dlmwrite('file_data_2g.csv',data2g,'-append');
%% Velocity simulation
%load('/media/turbik/DATA/PhD/MATLAB/Bacteria/c-stripes/data600.mat');
figure('Position',[600 500 620 350]);
load('/media/turbik/DATA/PhD/MATLAB/Bacteria/c-stripes/data_stationary_Dc=10.mat');
c = cp600(:,512)+cm600(:,512);
vxsim = (real(V_cm_600(:,512)).*cm600(:,512)+real(V_cp_600(:,512)).*cp600(:,512))./max(max(cp600(:,512)+cm600(:,512)));
vysim = (imag(V_cm_600(:,512)).*cm600(:,512)+imag(V_cp_600(:,512)).*cp600(:,512))./max(max(cp600(:,512)+cm600(:,512)));

vxsim = (real(V_cm_600).*cm600+real(V_cp_600).*cp600)/max(max(cp600+cm600))/max(max((real(V_cm_600).*cm600+real(V_cp_600).*cp600)/max(max(cp600+cm600))));
vysim = (imag(V_cm_600).*cm600+imag(V_cp_600).*cp600)/max(max(cp600+cm600))/max(max((real(V_cm_600).*cm600+real(V_cp_600).*cp600)/max(max(cp600+cm600))));

%c = c/max(c);

plot(linspace(0,960*2,1024*2),[vxsim(:,512);vxsim(:,512)],'LineWidth',2,'Color','b'); hold on;
plot(linspace(0,960*2,1024*2),[vysim(:,512);vysim(:,512)],'LineWidth',2,'Color','r');

data2h = [linspace(0,960*2,1024*2)',[vxsim(:,512);vxsim(:,512)],[vysim(:,512);vysim(:,512)]];

%axis tight;
xlabel(['{\ity} (',char(181),'m)']);
ylabel('{\itv_{x,y}}/{\itv}_0');
%legend({'{\itv^{+}_{x} + v^{-}_{x}}','{\itv^{+}_{x}c^+ + v^{-}_xc^-}'});
set(gca,'FontSize',18,'XTick',0:160:1000,'Units','pixels',...
    'Position',[125 100 420 250],'XMinorTick','on','YMinorTick','on');
%legend({'{\itv_{x} = v_x^+ c^+ + v_x^- c^-}','{\itv_{y}=v_y^+ c^+ + v_y^- c^-}'},'FontSize',14,'Box','off');
xlim([160+80 960-80]);
l = legend({'{\itv}_{x}','{\itv}_{y}'},'FontSize',14,'Box','off','Position',[0.7 0.69   0.1262    0.1125]);
set(gcf,'Position',[600 500 620 400],'PaperPositionMode','auto','Color','w');
ylim([-0.2 1.4]); %xlim([0 960+160]);
pause(3);
%print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/velocity_x,y_simulation','-dpng','-r300'); close(gcf);
%% data write into csv
fid = fopen('file_data_2h.csv','w'); 
fprintf(fid,'y(um),v_x(y)/v_0,v_y(y)/v_0\n');
fclose(fid);
dlmwrite('file_data_2h.csv',data2h,'-append');
%%
figure; plot(linspace(0,960,1024),imag(V_cp_600(:,512)),'b','LineWidth',2); hold on; 
plot(linspace(0,960,1024),imag(V_cm_600(:,512)),'r','LineWidth',2); 
plot(linspace(0,960,1024),cp600(:,512)*2.5,'.','Color','b','LineWidth',2);
plot(linspace(0,960,1024),cm600(:,512)*2.5,'.','Color','r','LineWidth',2);
legend({'{\itv_{y}^+}','{\itv_y^-}','2.5{\itc^+}','2.5{\itc^-}'},'FontSize',14,'Box','on');
set(gca,'FontSize',18,'XTick',0:160:1000);
xlabel(['{\ity}(',char(181),'m)']);
ylabel('{\itv_{x,y}/v}_0, {\itc^+, c^-}');
%% 
plot(linspace(0,960,1024),real(V_cm_600(:,1)).*cm600(:,1)+real(V_cp_600(:,1)).*cp600(:,1),'LineWidth',2); ylim([-0.25 0.5]); hold on;
plot(linspace(0,960,1024),imag(V_cm_600(:,1)).*cm600(:,1)+imag(V_cp_600(:,1)).*cp600(:,1),'LineWidth',2);
xlabel(['{\ity} (',char(181),'m)']);
ylabel('{\itv_y}/{\itv}_0');
legend({'{\itv^{+}_{x}c^+ + v^{-}_xc^-}','{\itv^{+}_{y}c^{+} + v^{-}_{y}c^{-}}'});
set(gca,'FontSize',14,'XTick',0:160:1000);

%% plotting v-map
sx = 10; sy = 10; ex = 10; ey = 10;
xs = 2;%2
ys = 5;%5
v_map = sqrt(vx_temp.^2+vy_temp.^2)*1.1906e-05/1.9725e-06;
image([X{1,1}(1,sx) X{1,1}(1,end-ex)]*1e6,[Y{1,1}(sy,1) Y{1,1}(end-ey,1)]*1e6,...
    (v_map(sx:end-ex,sy:end-ey)),'CDataMapping','scaled');
colormap(parula);
c = colorbar; 
c.Position = [0.9 0.3 0.03 0.26];
c.Limits = [0 25]/1e6;
c.Ticks = (0:5:20)*1e-6;
c.TickLabels = c.Ticks*1e6;
c.FontSize = 14;
c.Label.String = ['{\itv} (',char(181),'m/s)'];
%c.Label.FontSize = 14;
c.Label.Position = [0 32.5e-6];
c.Label.Rotation = 90;

%rectangle('Position',[20 900 300 100],'FaceColor','w','EdgeColor','w','LineWidth',2,'LineStyle','-');
%c.Ticks = (0:2:12);%c.TickLabels = {0:2:10};%c.LineWidth
hold on;
quiver(X{1,1}(sx:xs:end-ex,sy:ys:end-ey)*1e6,Y{1,1}(sx:xs:end-ex,sy:ys:end-ey)*1e6,...
    vx_temp(sx:xs:end-ex,sy:ys:end-ey),vy_temp(sx:xs:end-ex,sy:ys:end-ey),2,'w'); axis equal;
%axis equal; 
axis tight;

set(gca,'YDir','normal','FontSize',text_size,'Units','pixels',...
    'Position',[125 100 472 420],'YTick',(0:160:1000)-58,'YTickLabel',(0:160:960));
xlabel(['{\itx} (',char(181),'m)']);
ylabel(['{\ity} (',char(181),'m)']);

set(gcf,'Position',[600 600 620 520],'PaperPositionMode','auto','Color','w');

ylim([160+20 960-140]);
xlim([50 620]);
pause(2);

%% Cropped v-map
%figure('Position',[680 530 600 400],'PaperPositionMode','auto');
sx = 10; sy = 10; ex = 10; ey = 10;
xs = 2;%2
ys = 4;%5
v_map = sqrt(vx_temp.^2+vy_temp.^2)*1.1906e-05/1.9725e-06;
image([X{1,1}(1,sx) X{1,1}(1,end-ex)]*1e6,[Y{1,1}(sy,1) Y{1,1}(end-ey,1)]*1e6,...
    (v_map(sx:end-ex,sy:end-ey)),'CDataMapping','scaled');
%set(gcf,'Position',[680 530 560 500]);
colormap(parula);
hold on;
quiver(X{1,1}(sx:xs:end-ex,sy:ys:end-ey)*1e6,Y{1,1}(sx:xs:end-ex,sy:ys:end-ey)*1e6,...
    vx_temp(sx:xs:end-ex,sy:ys:end-ey),vy_temp(sx:xs:end-ex,sy:ys:end-ey),...
    1.5,'Color',[0.99 0.99 0.99],'LineWidth',2,'MaxHeadSize',2);
axis equal; axis tight;
set(gca,'YDir','normal','FontSize',text_size,'Units','pixels',...
    'Position',[0 0 420 630]);
axis off;
%xlabel(['{\itx} (',char(181),'m)'],'FontSize',text_size);
%ylabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
ylim([200 500]);
xlim([380 580]);

set(gcf,'Position',[600 500 420 630],'PaperPositionMode','auto');
print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/velocity_cropped_map','-dpng','-r600');
%% Histogram of velocity
figure('Position',[680 530 600 400],'PaperPositionMode','auto');
vx_temp1 = vx_temp*1.1906e-05/1.9725e-06;
histogram(vx_temp1,linspace(-2e-5,2e-5,50),'FaceColor','b','DisplayStyle','stairs','LineWidth',2,'Normalization','probability'); 
hold on;
histogram(vy_temp*1.1906e-05/1.9725e-06,linspace(-2e-5,2e-5,50),'FaceColor','r','DisplayStyle','stairs','LineWidth',2,'Normalization','probability');
%set(gcf,'Position',[680 530 560 420],'PaperPositionMode','auto');
%axis square;
legend({'x','y'},'Box','off');
set(gca,'FontSize',18,'YScale', 'log','XLim',[-10 20]*1e-6,'XTick',(-15:5:20)*1e-6,'XTickLabel',(-15:5:20),'TickLength',[0.02 0.5],...
    'Unit','pixel','Position',[125 100 450 250],'XMinorTick','on');
%ylabel(['{\ity} (',char(181),'m)'],'FontSize',16);
xlabel(['\langle{\itv}\rangle_{\itt} (',char(181),'m/s)'],'FontSize',18);
ylabel('Probability','FontSize',18); box('on');
ylim([1e-3 1]);
%% Misha velocity
figure('Position',[680 530 600 400],'PaperPositionMode','auto');
%linspace(0, 960, 1024);
%plot(linspace(0, 960, 1024)',real(V_cm_600(:,1)),'LineWidth',2,'Color','b'); hold on; plot(linspace(0, 960, 1024)',real(V_cp_600(:,1)),'LineWidth',2,'Color','b','LineStyle','--');
%plot(linspace(0, 960, 1024)',imag(V_cm_600(:,1)),'LineWidth',2,'Color','r'); hold on; plot(linspace(0, 960, 1024)',imag(V_cp_600(:,1)),'LineWidth',2,'Color','r','LineStyle','--');
plot([linspace(0, 960, 1024)'],[real(V_cm_600(:,1))+real(V_cp_600(:,1))],'LineWidth',2,'Color','b'); hold on;
plot([linspace(0, 960, 1024)'],[imag(V_cm_600(:,1))+imag(V_cp_600(:,1))],'LineWidth',2,'Color','r');
%plot([linspace(0, 960, 1024)',linspace(0, 960, 1024)'],[imag(V_cm_600(:,1)),imag(V_cp_600(:,1))],'LineWidth',2);
%plot([imag(V_cm_600(:,1)),imag(V_cp_600(:,1))]);

text_size = 18;

%plot(fliplr(yax(1:end-astart)),vx_y(1:end-astart)/1e-6,'.-','LineWidth',1.5);
%hold on; plot(fliplr(yax(1:end-astart)),vy_y(1:end-astart)/1e-6,'.-','LineWidth',1.5);
%axis square;

axis tight;
set(gca,'FontSize',text_size,'YLim',[-0.25 0.5],'XTick',0:160:1000,...
    'Unit','pixel','Position',[125 100 420 250],'XMinorTick','on','YMinorTick','on'); 
legend({'{\itc}^{\fontsize{8}+}','{\itc}^-'},'Box','off','Location','best');
xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
%ylabel(['{\itv}_{{\itx},{\ity}}({\ity})/{\itv}_0'],'FontSize',text_size);
%ylabel(['{\itv}_{\itx} / {\itv}_0'],'FontSize',text_size);
%ylabel(['{\itv}_{\ity} / {\itv}_0'],'FontSize',text_size);
ylabel(['{\itv}/{\itv}_0'],'FontSize',text_size);
legend({'{\itv}_{\itx}','{\itv}_{\ity}'},'Box','off','Location','best');
%print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/sim_velocity_x,y_2','-dpng','-r300');
%% Trajectories from TrackMate Fiji
filename = '/media/turbik/DATA/PhD/DATA/2017_Bacteria_c,s,v-stripes/NewRM_cells/Cell_3_5_for_analysis/All Spots statistics.csv';
fileID = fopen(filename,'r');

dataArray = textscan(fileID, '%*s%*s%s%*s%s%s%*s%*s%s%[^\n\r]', 'Delimiter', ',',  'ReturnOnError', false);

fclose(fileID);

raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN};
data = raw;
%%
B = cell2mat(data);
[~,idx] = sort(B(:,1));
C = B(idx,:);
ID = C(:,1);
X = C(:,2);
Y = C(:,3);

maC = max(ID);

region = cell(maC,3);
for i = 1:maC
    region{i,1} = find(C(:,1) == i);
    region{i,2} = X(region{i,1});
    region{i,3} = Y(region{i,1});
    
    plot(region{i,2},region{i,3},'.'); hold on;
    axis equal;
    
end
plot(X(ID == 847),Y(ID == 847));
%% plotting
longTraj = cell(1);
k = 1;
for i = 60:maC
    if max(region{i,3})-min(region{i,3}) > 400
        plot(region{i,2},region{i,3},'.'); hold on;
        pause(0.2); title(num2str(i));
        longTraj{k} = num2str(i);
        k = k + 1;
        legend(longTraj);
    end
    axis equal;
end
%% velocity vs curvature of the pattern
L = 160; %µm
scale = 8192/262.02/10*0.97; 
% longTraj 

%k = longTraj; %, 60, 65, 107, 847, 163
y0 = 610;%1096, 596, 126

for i = 1:length(longTraj)-4
    
s = str2num(longTraj{i});
x = region{s,2};
y = region{s,3};

if s == 65
    x = region{s,2}(200:end);
    y = region{s,3}(200:end);
elseif s == 107
    x = region{s,2}(1:end-100);
    y = region{s,3}(1:end-100);
elseif s == 163
    x = region{s,2}(1:end-220);
    y = region{s,3}(1:end-220);
end
y = y0 - y;

x = smooth(x/scale,20);
y = smooth(y/scale,20);

col = rand(1,3);

subplot(1,2,1);
hold on;
[X,Y] = meshgrid(0:10:400, -380:10:180); 
quiver(X,Y,cos(pi*Y/160),sin(-pi*Y/160),0.25,'Color',[0,113,188]/256,'ShowArrowHead','off');
quiver(X,Y,-cos(pi*Y/160),-sin(-pi*Y/160),0.25,'Color',[0,113,188]/256,'ShowArrowHead','off');
plot(x,y,'.','Color',col);

nx = cos(pi*y(2:end)/L);
ny = sin(pi*y(2:end)/L);

vx = diff(x); 
vy = diff(y);
xlabel(['{\itx} (',char(181),'m)'],'FontSize',text_size);
ylabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
quiver(x(2:2:end), y(2:2:end), vx(1:2:end), vy(1:2:end),2,'ShowArrowHead','on','Color',col);%,'Color','r'
axis equal;
set(gca, 'FontSize',text_size);

subplot(1,2,2);
deltaTh = atan(vx./vy)-(-atan(nx./ny));
hold on;
deltaTh(abs(deltaTh) > pi/2) = deltaTh(abs(deltaTh) > pi/2) - pi;
plot((y(1:end-1,1)), smooth(deltaTh*180/pi),'.','Color',col);
%plot(-atan((cos(pi*smooth(y(3:end),20)/scale/160)./sin(pi*smooth(y(3:end),20)/scale/160))),'.');
ylim([-90 90]);
%set(gca);
set(gca, 'FontSize',text_size);

xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
ylabel(['{\it\theta_v - \theta_n} (^o)'],'FontSize',text_size);
end

%% Bend instab 10um thickness
figure('Position',[680 530 1392 274]); movegui(gcf,'center');
plot(H/(1392/688.17),D/(160e-6*10e-6*688.17e-6),'o-'); xlabel(['x (',char(181),'m)']); ylabel('c (m^{-3})');
%%
[~,nn] = size(results);
pfft = cell(nn,1);
mfft = zeros(nn,1);
for i = 1:nn
    clear A;
    A = nanmean(results{4,i});
    A=(A(~isnan(A)));
    pfft{i} = A;
    gf = fit((1:length(pfft{i}))',pfft{i}','fourier1');
    mfft(i) = gf.w;
end
%%
for i = 1:nn
    A = nanmean(results{4,i});
    hold on;
    plot3(ones(length(A),1)*i/5,1:length(A),A);
    view(-96,24);
    grid on; 
    xlabel('t(s)'); ylabel(['x (',char(181),'m)']); zlabel(['v_{y} (',char(181),'m/s)']);
end
%%
con = zeros(921,1);
for k =1:921
    %s = dir(['/media/turbik/DATA/PhD/DATA/2018_Bacteria_cicr_+1/results_coordinates/Results_',num2str(k),'.csv']);
    %if s.bytes < 1000
    %    D(:,k) = nan;
    %    continue
    %end
    data = dlmread(['Cell_fluoresc_5fps/Results/',...
        'Results_',num2str(k),'.csv'],',',1,2);    
    con(k) = length(data(:,7));
end
%% 3D plot image
Time = [0,nFiles]/5;
Radius = 150-([H(end) H(1)])/scale;
figure('Position',[680 530 600 420]);
image(Time,Radius,D/dV*1e-14,'CDataMapping','scaled'); 
%num2str((sort(get(gca,'YTick'),'descend'))','%.0f')
set(gca, 'YTickLabel',{150,100,50,0},'FontSize',text_size,...
    'Units','pixels', 'Position', [87.5667 67.3000 429.2833 349.9500]);

c = colorbar; c.Label.String = '{\itc} (\times10^{14} m^{-3})';
%c.Label.FontSize = 14;
%c.Ticks = (0:max(max(D/dV))*1e-14);
%c.Label.String = ['Veclocity (',char(181),'m/s)'];
c.Label.FontSize = text_size;
%c.LineWidth = 2;
%c.Label.Position = [-0.5 -3];
c.Label.Position = [-0.5 -5];
%c.Position = [0.84 0.475 0.03 0.4];
c.Position = [0.88 0.5 0.03 0.4];
%c.Ticks = (0:2:12);
%c.TickLabels = {0:2:10};
%c.Label.Rotation = 0;
%c.LineWidth

%view(0,90); 
%set(gca,'FontSize',text_size,'LineWidth',1); grid off; 
xlabel('Time (sec)','FontSize',text_size); %xlim([0 50])
ylabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
%box on; axis square; 
%xlim([0 max(max(Time))]); 
%ylim([0 max(max(Radius))]); 
%ylim([0 1000]); %axis tight;
%title('Bacteria number distribution (lower conc)');
%colorbar;%('Ticks',[0,1,2,3,4]);
%% Average velocity (for the case of 30 step)
tempx = nanmean(results{3,101},2)/scale/0.1*1e6;
tempy = (max(results{2,1}(:,1))-results{2,1}(:,1))/scale*1e6;
plot(-tempx(5:end-10,1),tempy(5:end-10,1),'o');
axis square;
xlabel(['v_{x} (',char(181),'m/s)'],'FontSize',text_size); 
ylabel(['y (',char(181),'m)'],'FontSize',text_size);
set(gca,'FontSize',text_size,'XLim',[0 8],'YLim',[0 1000]);
%,nanmean(results{4,601},2)/scale/0.1*1e-6,results{2,1}(:,1)/scale*1e-6,'ro');
%%
res = cell(2,1);
res{1} = zeros(311,231);
res{2} = zeros(311,231);
for h = 1:311
    for w = 1:231
        avg1 = 0;
        avg2 = 0;
        for i = 1:1199
            if ~isnan(results{3,i}(h,w))
                res{1}(h,w) = res{1}(h,w) + results{3,i}(h,w);
                avg1 = avg1 + 1;
            end
            if ~isnan(results{4,i}(h,w))
                res{2}(h,w) = res{2}(h,w) + results{4,i}(h,w);
                avg2 = avg2 + 1;
            end
        end
        res{1}(h,w) = res{1}(h,w)/avg1;
        res{2}(h,w) = res{2}(h,w)/avg2;
    end
end

%% Histograms of velocity components
hold on;
histogram(-(results{3,101}/scale/0.2),(-1e-5:0.5e-6:2e-5)/1e-6,'FaceColor','b','DisplayStyle','stairs','LineWidth',2); xlim([-1.5e-5 1.5e-5]);axis square;
histogram(-(results{4,101}/scale/0.2),(-1e-5:0.5e-6:2e-5)/1e-6,'FaceColor','r','DisplayStyle','stairs','LineWidth',2); xlim([-1.5e-5 1.5e-5]);axis square;
%subplot(1,2,2); 
%axis square;
axis tight;
legend({'x','y'},'Box','off');
set(gca,'FontSize',18,'XLim',[-12,12],'XTick',-10:5:10,'TickLength',[0.02 0.5]);
%ylabel(['{\ity} (',char(181),'m)'],'FontSize',16);
xlabel(['{\itv} (',char(181),'m/s)'],'FontSize',18);
ylabel('Count','FontSize',18); box('on');
%%
hold on;
histogram(-(results{3,101}(:,170:180)/scale/0.2),(-1e-5:0.5e-6:2e-5)/1e-6,'FaceColor','b','DisplayStyle','stairs','LineWidth',2); xlim([-1.5e-5 1.5e-5]);axis square;
histogram(-(results{4,101}(:,170:180)/scale/0.2),(-1e-5:0.5e-6:2e-5)/1e-6,'FaceColor','r','DisplayStyle','stairs','LineWidth',2); xlim([-1.5e-5 1.5e-5]);axis square;
%subplot(1,2,2); 
%'Normalization','probability',
%axis square;
axis tight;
legend({'x','y'},'Box','off')
set(gca,'FontSize',18,'XLim',[-12,12],'XTick',-10:5:10,'TickLength',[0.02 0.5]);
%ylabel(['{\ity} (',char(181),'m)'],'FontSize',16);
xlabel(['{\itv} (',char(181),'m/s)'],'FontSize',18);
ylabel('Count','FontSize',18); box('on');
%% Vectors
scale = 8192/262.02/10;
text_size = 18;
g=10; a = 20; b = 19;
max_y = max(max(results{2,601}));
max_vy = max(max(results{4,601}(a:g:end-b,a:g:end-b)));
figure('Position',[200 1000 1200 800]);
quiver(results{1,601}(a:g:end-b,a:g:end-b)/scale,...
    (max_y-results{2,601}(a:g:end-b,a:g:end-b))/scale,...
    results{3,601}(a:g:end-b,a:g:end-b),...
    (-results{4,601}(a:g:end-b,a:g:end-b)),2,'Color',[0 0.7 0]);
%axis([0 750 0 1500]);
%g=5;image([0,0.75e-3],[0,1e-3],sqrt(results{3,601}.^2+results{4,601}.^2)/max(max(results{3,601}))*256); axis([0 0.75e-3 0 1e-3]); colormap(gray)
%set(gca,TickLabelFormat,'$%,.0f')
axis equal;
axis tight;
set(gca,'FontSize',text_size/2,'XLim',[0 750],'YLim',[0 1000],...
    'XTick',0:200:800,'YTick',0:200:1000);
ylabel(['y (',char(181),'m)'],'FontSize',text_size/2);
xlabel(['x (',char(181),'m)'],'FontSize',text_size/2);
box('on');

%%
k = 101;
g=1; a1 = 8; b1 = 18;
gim = 3; a = 8; b = 18;

figure;
imh = image([min(results{1,k}(1,a:gim:end-b)) max(results{1,k}(1,a:gim:end-b))]/scale,...
    [min(results{2,k}(a:gim:end-b,1)) max(results{2,k}(a:gim:end-b,1))]/scale,fliplr(sqrt((results{3,k}(a1:g:end-b1,a1:g:end-b1)).^2+...
    (results{4,k}(a1:g:end-b1,a1:g:end-b1)).^2)/max(max(sqrt((results{3,k}(a1:g:end-b1,a1:g:end-b1)).^2+(results{4,k}(a1:g:end-b1,a1:g:end-b1)).^2)))*64));
hold on;
pos = get(gca,'Position');
pos(2) = pos(2) - 0.01;

quiver((2100-results{1,k}(a:gim:end-b,a:gim:end-b))/scale,...
    (results{2,k}(a:gim:end-b,a:gim:end-b))/scale,...
    -results{3,k}(a:gim:end-b,a:gim:end-b),...
    (results{4,k}(a:gim:end-b,a:gim:end-b)),2,'Color','w');
set(gca,'Position',pos);
colormap(parula);
c = colorbar('northoutside');
c.Label.String = ['{\itv} (',char(181),'m/s)'];
c.Label.FontSize = text_size/2;
c.Label.Position = [-5 -0.1];
c.Ticks = linspace(0,64,6);
c.TickLabels = {0:2:10};
c.Position = [0.4 0.96 0.3 0.03];
c.XAxisLocation = 'bottom';
axis equal;
axis tight;
set(gca,'FontSize',text_size/2,'XLim',[20 650],'YLim',[20 900],...
    'XTick',20:200:620,'XTickLabel',0:200:800,'YTick',20:200:1220,'YTickLabel',0:200:800,'YDir','normal');
ylabel(['{\ity} (',char(181),'m)'],'FontSize',text_size/2);
xlabel(['{\itx} (',char(181),'m)'],'FontSize',text_size/2);
rectangle('Position',[(2300-1700)/scale 0 10/scale 3100/scale],'EdgeColor','w');

%% vx
d = nanmean(results{3,k}(:,:),2)/scale/0.2;
%d(isnan(d)) = 0;
plot(results{2,k}(:,1)/scale,-d,'o-');
%,'XLim',[20 650],'YLim',[20 900],...
    %'XTick',20:200:620,'XTickLabel',0:200:800,'YTick',20:200:1220,'YTickLabel',0:200:800,'YDir','normal');
%ylim([-5 10]);
ylabel(['{\itv_{x}} (',char(181),'m/s)'],'FontSize',text_size);
xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
set(gca,'FontSize',text_size,'XTick',0:200:1000);
print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/veloc_dist_vx_tot','-dpng','-r0');

%% vy
d = nanmean(results{4,k}(:,170:180),2)/scale/0.2;
%d(isnan(d)) = 0;
plot(results{2,k}(:,1)/scale,d,'o-');
%,'XLim',[20 650],'YLim',[20 900],...
    %'XTick',20:200:620,'XTickLabel',0:200:800,'YTick',20:200:1220,'YTickLabel',0:200:800,'YDir','normal');
ylim([-5 5]);
ylabel(['{\itv_{y}} (',char(181),'m/s)'],'FontSize',text_size);
xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
set(gca,'FontSize',text_size,'XTick',0:200:1000);
print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/veloc_dist_vy','-dpng','-r0');
%% |v|
d = sqrt((nanmean(results{3,k}(:,170:180),2)/scale/0.2).^2 + (nanmean(results{4,k}(:,170:180),2)/scale/0.2).^2);
%d(isnan(d)) = 0;
plot(results{2,k}(:,1)/scale,d,'o-');
%,'XLim',[20 650],'YLim',[20 900],...
    %'XTick',20:200:620,'XTickLabel',0:200:800,'YTick',20:200:1220,'YTickLabel',0:200:800,'YDir','normal');
ylabel(['|{\itv}| (',char(181),'m/s)'],'FontSize',text_size);
xlabel(['{\ity} (',char(181),'m)'],'FontSize',text_size);
set(gca,'FontSize',text_size,'XTick',0:200:1000);
print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/veloc_dist_absv','-dpng','-r0');


%% Traj continue
%x = smooth(x,10); y = y;
for i=1:length(x)-2
    vx(i) = (x(i+1) - x(i));
    vy(i) = (y(i+1) - y(i));
end
plot(y(1:end-2),vx); hold on; plot(y(1:end-2),vy); xlim([0 1000])
%% concentration fit
mn = mean(D,2)/dV;
Y = mn(1000:2000,1)*1e12;
X =((H(1000:2000))/scale/1e6)';%(-90+(1:length(Y))');
%Y = mn(36:50,1); X = 600-(1000-H(36:50)/scale)';
plot(X(1:end),Y(1:end)/1e14,'o','MarkerFaceColor','b');
hold on;
%a = 1.933e-10; b = 74.27; c = 22.84; x0 = -0.131; area > 60;
%a = 1.933e-10; b = 74.57; c = 22.84/scale*1e-6; x0 = 2.131;
a = 5e0; b = 0.2; L = 150/1e6; x0 = -0.1246;
%a1 =   5.002e+14/1e14; b1 = 0.3499; c1 = 4.299;
%x = X;
%x = (-100:0.1:250);
%plot(x,a*exp(b*L/pi*cos(pi*x/L))+a*exp(-b*L/pi*cos(pi*x/L)));
%plot(x,a1*exp(-((x-b1)/c1).^2),'LineWidth',2);%cos(pi*(x-x0)/L)
%axis square;
%set(gca,'FontSize',18,'XLim',[-50,50],'XTick',-50:25:50,'TickLength',[0.02 0.5]);
%xlim([-50 50]); legend('Data','Gaussian fit');
%xlabel(['y (',char(181),'m)'],'FontSize',18);
%ylabel('c (\times 10^{14} m^{-3})','FontSize',18); box('on');
%% Bend instab
t = (0:921)/5;
scale = 8192/262.02/10; %px/um --> 20x (DF, 2.1mm WD)
h = 508; w = 2100;
nFiles = 922;%922
totConc = zeros(nFiles,1);
test = zeros(nFiles,2);
dStep = 5*scale;
H = 0:dStep:h;
nDist = length(H);
D = zeros(nDist,nFiles);
%normD = zeros(nDist,nFiles);
%meanD = zeros(nFiles,1);
%
dstp = floor(linspace(10,50,922));
for k = 1:922
    %if s.bytes < 1000
    %    D(:,k) = nan;
    %    continue
    %end
    data = dlmread(['Cell_fluoresc_5fps/Results/',...
        'Results_',num2str(k),'.csv'],',',1,2); 
    %'Cell_3_5_for_analysis/Results_',num2str(k),'.csv'],',',1,2);
    X = data(:,7); Y = data(:,8); area = data(:,1); AR = data(:,32); theta = data(:,18);
    cond = Y>150 & area > 10 & AR > 1;
    x = X(cond); y = Y(cond); theta = theta(cond);
    theta(theta>90) = theta(theta>90) - 180;
    N(k) = length(x);
    
    dx = linspace(min(x),max(x),dstp(k));
    dw1 = 0; dw2 = 0; ddw = 0; dtheta = 0; 
    for j = 1:length(dx)-1
        cond = x>=dx(j) & x<=dx(j+1);
        if length(nonzeros(cond)) == 0
            dw1(j) = nan;
            dw2(j) = nan;
            dtheta(j) = nan;
        else
            dw1(j) = max(y(cond));
            dw2(j) = min(y(cond));
            dtheta(j) = mean(theta(cond));
        end
    end
    
%     for j = find(isnan(dw1))
%         dw1(j) = mean([dw1(j-1),dw1(j+1)]);
%         dw2(j) = mean([dw2(j-1),dw2(j+1)]);
%         dtheta(j) = mean([dtheta(j-1),dtheta(j+1)]);
%     end
    
    dw1(isnan(dw1)) = 0;
    dw2(isnan(dw2)) = 0;
    dtheta(isnan(dtheta)) = 0;
    
    ddw = nanmean(dw1 - dw2)*cosd(dtheta);    
    
    W(k) = nanmean(ddw);
    %L(k) = sum(diff(dx))*cosd(nanmean(dtheta));
    %L(k) = sum(diff(dx))/cosd(std(theta)*exp(1)/2);
    
    %nansum(sqrt(1+tand(dtheta).^2).*diff(dx))
    %L(k) = nansum(sqrt(1+tan(dtheta*pi/180).^2).*diff(dx));
    %L(k) = nansum(diff(dx)./(cos(dtheta*pi/180)));
    L(k) = nansum(diff(dx)'./smooth(cos(dtheta*pi/180)));
    disp(L(k));
    dV1(k) = mean((diff(dx)'./smooth(cos(dtheta*pi/180)))'.*(dw1 - dw2).*cosd(dtheta))*20e-18/2;
    %L(k) = mean(sqrt(1+tand(dtheta).^2)+sum(diff(dx)));
    %w/cosd(std(theta));
%     
    %pause(2);
    test(k,1) = std(y);
    theta(theta>90) = theta(theta>90)-180;
    test(k,2) = std(theta);
    
    T(k) = std(theta);
    
    test(k,3) = max(theta);
%   if k < 922
%       test(k,2) = nanmean(nanmean(sqrt(results{4,k}.^2+results{3,k}.^2)));
%       test(k,3) = nanmean(nanmean(results{3,k}));
%       test(k,4) = nanmean(nanmean(results{4,k}));
%   end
    area = data(:,1);
    nPoints = length(X);
    for i = 1:nPoints
        if X(i)>0 && area(i)>0 && Y(i)>135
            totConc(k) = totConc(k) + 1;
            for j = 1:nDist-1
                if Y(i) > H(j) && Y(i) < H(j+1)
                    D(j,k) = D(j,k) + 1;
                end
            end
        end
    end
end
delta_y = h/nDist;
dV = delta_y*w/scale^2*20*1e-18;
text_size = 18;

%plot(N); hold on; plot(W); plot(L-w,'.');
%legend({'Number of bacteria','Width','Length'},'FontSize',12,'Location','best');
%% plOTTTING UNDULATION FITTING
im = imread(['Cell_fluoresc_5fps/Result of FormatFactory20x_cell2_',...
    '5fps_2-10921.jpg']);
image(im); colormap(gray(256)); hold on; plot(250-50*sin((0:2000)*pi/160),'r-'); axis equal;
f = fittype('a*sin(x*pi/160)');
fit1 = fit(x,y,f,'StartPoint',35,'Upper',[]);
%plot(250-fit1(1:2000),'b');
df = differentiate(fit1,1:2000);
sum(sqrt(1+df.^2));

%% plotting c map

i = 3;
image(cp{i,1}+cm{i,1},'CDataMapping','scale'); colormap(parula);
c = colorbar;
set(gca,'Units','pixels','Position',[125 100 420 250],'FontSize',18);
set(gcf,'Position',[500,500,620,350],'PaperPositionMode','auto','Color','w');

c.Location = 'south';
c.Position = [0.5 0.2 0.3 0.05];
c.FontSize = 14;
c.Limits = [0 1.2];%max(max(cp{i,1}+cm{i,1}))];
c.Ticks =0:0.4:1.6;
%c.Ticks =0:0.4:1.6;
%c.TickLabels = 12;
c.Label.String = '{\itc}/\langle{\itc}\rangle';
c.Label.Position = [-0.2 1.25];

hold on; axis equal; axis off;
inval = 1:20:1024;
[x,y] = meshgrid(inval,inval);
quiver(x,y,cos(theta{i,1}(inval,inval)),sin(theta{i,1}(inval,inval)),0.25,'Color',[1,0.25,0.25],'ShowArrowHead','off','LineWidth',1.5);
quiver(x,y,-cos(theta{i,1}(inval,inval)),-sin(theta{i,1}(inval,inval)),0.25,'Color',[1,0.25,0.25],'ShowArrowHead','off','LineWidth',1.5);
xlim([0 768]); ylim([512-256/1.5 512+256/1.5]);

text(25,625,['{\itt/T} = ',num2str(k(i)/10000)],'Color','k','BackgroundColor','w','FontSize',14);

pause(2);

%% max concentration
%k = [50:50:500,600:100:900,1000:500:10000];
k = [50:50:10000];
l = length(k);
for i = 1:l
    [~,~,~,var4,var5]= ReadData2(foldname,k(i),1024);
    i
    cmax(i) = max(max((var4+var5)));
end

%% plotting
for i = 1:4
    figure;
    image(cp{i,1}+cm{i,1},'CDataMapping','scale')
end

%% plotting cmax
plot(k/max(k),cmax/0.08,'bo-','MarkerFaceColor','b'); hold on;
%plot([0.4333,0.4333],[0,20],'r--','LineWidth',2);

data3i = [k'/max(k),cmax'/0.08];
xlabel('{\itt/T}'); ylabel('{\itc}_{max}/\langle{\itc}\rangle');
ylim([0 25]);

set(gca,'Units','pixels','Position',[125 100 360 250],'FontSize',18,'XTick',0:0.2:1);
set(gcf,'Position',[500,500,620,350],'PaperPositionMode','auto','Color','w');
pause(2);
%% plotting c spreading
cspread = zeros(1024,1024+4);
for i = 1:1024
    cspread(:,i+4) = (cp{4,1}(:,i)+cm{4,1}(:,i))/0.08;
end
for i = 1:3
    cspread(:,i) = (cp{i,1}(:,512)+cm{i,1}(:,512))/0.08;
end
figure('Position',[500,500,620,380],'PaperPositionMode','auto','Color','w');
x = linspace(0,320,1024);
h1 = plot(x,cspread(:,4:end),'b','LineWidth',2); hold on;
h4 = plot(x,cspread(:,3),'k','LineWidth',2);
h3 = plot(x,cspread(:,2),'r','LineWidth',2);
h2 = plot(x,cspread(:,1),'Color',[0 0.75 0],'LineWidth',2);

data3j = [x',cspread(:,1),cspread(:,2),cspread(:,3),cspread(:,5:64:end)];

%pause(2);
%plot(k/max(k),cmax/0.08,'bo-'); hold on;
%plot([0.4333,0.4333],[0,20],'r--','LineWidth',2);

xlabel(['{\ity} (',char(181),'m)']); ylabel('{\itc}/\langle{\itc}\rangle');
xlim([120 200]); ylim([0 20]);
legend([h2,h3,h4,h1(1)],{'0','0.03','0.1','1'},'FontSize',14,'Box','off');
%title(hleg,'{\itt/T}');
%hleg.Title.Visible = 'on';
%hleg.Title.NodeChildren.Position = [0.5 1.5 0];
set(gca,'Units','pixels','Position',[125 100 280 250],'FontSize',18,'XTick',80:20:320);

pause(2);
%close(gcf);

%% data into csv
fid = fopen('file_data_3i.csv','w'); 
fprintf(fid,'t/T,c_max/<c>\n'); 
fclose(fid);
dlmwrite('file_data_3i.csv',data3i,'-append');

str0 = 'c/<c>(t/T=1;';
str1 = '';
for jjj = 64:64:1024
    str1 = [str1,str0,'x=',num2str(jjj/1024),'*x_max),'];
end
str1 = str1(1:end-1);

fid = fopen('file_data_3j.csv','w'); 
fprintf(fid,['y(um),c/<c>(t/T=0),c/<c>(t/T=0.03),c/<c>(t/T=0.1),',str1,'\n']);
fclose(fid);
dlmwrite('file_data_3j.csv',data3j,'-append');
%%

figure('Position',[600 500 620 350],'PaperPositionMode','auto','Color','w');
cp60=cp60(:,512); cm60=cm60(:,512); cp600=cp600(:,512); cm600=cm600(:,512);
h1 = plot(linspace(1,320,1024),c(:,:),'b','LineWidth',1); hold on;

h2 = plot(linspace(1,960,1024),(cp600+cm600)*2,'r','LineWidth',2);
%plot([1,1],[1,2],'r','LineWidth',2); hold on;
%plot([1,1],[1,2],'b','LineWidth',2);

ylabel('{\itc}/{\itc}_{0}');
xlabel(['{\ity} (',char(181),'m)']);

xlim([80 160+80]); %ylim([0 0.6])
set(gca,'Units','pixels','Position',[125 100 420 250],'FontSize',18,'XTick',80:40:240); 

legend([h2,h1(1)],{'Straight jet','Undulating jet'},'FontSize',14);
%set(gcf,'Position',[675,524,570,350],'PaperPositionMode','auto');

% plot(H/scale+76,D1(:,[10,300,600,910]),'LineWidth',2);

plot([150 170],[0.4 0.3],'r','LineWidth',2); hold on;
plot([150 170],[0.4 0.3],'b','LineWidth',2);

legend({'2 s','60 s','120 s','180 s'}); 
xlabel(['{\ity} (',char(181),'m)']); ylabel('Normalized concentration');
%ylabel('{\itc} (\times10^{14} m^{-3})');
set(gca,'FontSize',14,'XTick',80:20:240); xlim([100 160+60]); ylim([0 0.6])
set(gcf,'Position',[375,524,570,350],'PaperPositionMode','auto');
legend({'before wave','during wave'});

%subplot(1,2,2);
%% colorbar 0 0.6
figure('Color','w');
image([0 0.5; 0 0]); colormap(parula(300)); c = colorbar;
set(gca,'Position',[0.1 0.1 0.1 0.1]);
c.Limits = [0 300]; c.Location = 'north';
c.TickLabels = {'0','0.2','0.4','0.6'};%,'0.4','0.5','0.6'
c.FontSize = 14; c.Units = 'pixel'; c.Position = [100 200 200 20];
%%
x11 = (1:922)/5;
y11 = smooth(max(D1),5)./cosd(test(:,2));
%y11 = sum(D)'/(20e-6*20e-6*w*1e-6);
plot(x11(1:20:end),y11(1:20:end),'bo','LineWidth',1,'MarkerFaceColor','b');
xlabel('Time (s)'); ylabel('{\itc} (\times10^{14} m^{-3})'); 
xlim([0 180]);
%set(gcf,'Color','white','Position',[300 600 1200 350]);
set(gcf,'Position',[375,524,570,350],'PaperPositionMode','auto');
%set(gcf,'Position',[375,524,400,350],'PaperPositionMode','auto');
set(gca,'FontSize',14,'XTick',0:25:200);
%% Number vs time
fr = sum(D,1);
plot((1:20:size(D,2))/5, fr(1,1:20:end),'bo','LineWidth',1);
set(gca,'FontSize',14,'XLim',[0 180],'YLim',[50 200]);
ylabel('Count');
xlabel('time (s)');
set(gcf,'Position',[675,524,570,350],'PaperPositionMode','auto');
%print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/number_undulat','-dpng','-r300');

%% Velocity
scale = 8192/262.02/10; %px/um --> 20x (DF, 2.1mm WD)
h = 3100; w = 2300;
nFiles = 921;
totConc = zeros(nFiles,1);
test = zeros(nFiles,2);
nDist = floor(h/10);
H = linspace(0,h,nDist);
%length(dir([folder,'*.jpg']));
%D = zeros(nDist,nFiles);
normD = zeros(nDist,nFiles);
meanD = zeros(nFiles,1);

for k = 1:921
        test(k,2) = nanmean(nanmean(sqrt(results{4,k}.^2+results{3,k}.^2)));
        test(k,3) = nanmean(nanmean(results{3,k}));
        test(k,4) = nanmean(nanmean(results{4,k}));
end

delta_y = h/nDist;
dV = delta_y*w/scale^2 * 20 * 1e-18;
text_size = 18;
plot((1:10:920)'/5,(test(1:10:end-1,2)/scale/0.2),'ro','LineWidth',1); hold on;
%plot(t,T*pi/180);
ylim([0 20]); xlim([0 180]);
xlabel('{\itt} (s)','FontSize',18);
ylabel(['{\itv} (',char(181),'m/s)'],'FontSize',18);
set(gca,'FontSize',18,'Units','pixels','Position',[125 100 420 300],'Box','on');
set(gcf,'Position',[600,500,620,400],'PaperPositionMode','auto');
%% Number bacteria & Velocity on one plot
%plot((1:20:size(D,2))/5, fr(1,1:20:end),'bo','LineWidth',1); hold on;
%plot((1:10:920)'/5,(test(1:10:end-1,2)/scale/0.2),'ro','LineWidth',1);

d1 = 1;%20
d2 = 1;%10
t1 = (1:d1:size(D,2))/5;
t2 = (1:d2:920)'/5;

y1 = smooth(fr(1,1:d1:end),5);
y2 = (test(1:d2:end-1,2)/scale/0.2);

[ax,h1,h2] = plotyy(t1(1:20:end),y1(1:20:end),t2,y2);
set(h1,'LineStyle','none','Marker','o','Color','b','MarkerFaceColor','b');
set(h2,'LineStyle','none','Marker','o','Color','r','MarkerFaceColor','r');

set(ax(1),'FontSize',14,'XLim',[0 180],'YLim',[50 200]);
set(ax(2),'FontSize',14,'XLim',[0 180],'YLim',[0 20]);
%%
%plot(((1:length(totConc))-1)'/5,totConc,'.');
figure;
plot(((1:length(totConc))-1)'/5,test(:,1)/scale,'.r');
%plot(((1:length(totConc))-1)'/5,test(:,2)/scale/0.2,'.g');
axis tight;
%legend({'x','y'},'Box','off')
set(gca,'FontSize',18)%,'XLim',[-10,10],'XTick',-10:5:10,'TickLength',[0.02 0.5]);
%ylabel(['{\ity} (',char(181),'m)'],'FontSize',16);
xlabel('Time (s)','FontSize',18);
%ylabel('Number of bacteria','FontSize',18); box('on');
ylabel(['Width of distortion ',' (',char(181),'m)'],'FontSize',18); box('on');
%ylabel(['Velocity magnitude',' (',char(181),'m/s)'],'FontSize',18); box('on');
%%
text_size = 18;
figure('Position',[680 530 610 420]);
[ax,pl1,pl2] = plotyy(((1:length(totConc))-1)'/5,test(:,1)/scale,((1:length(totConc)-1)-1)'/5,smooth(test(1:end-1,2)/scale/0.2,100));
%axis tight;
xlabel('Time (s)','FontSize',text_size);
ylabel(ax(1),['max(2{\itW}_{0}) ',' (',char(181),'m)'],'FontSize',text_size,'Color','r');
ylabel(ax(2),['\langle{\itv}\rangle',' (',char(181),'m/s)'],'FontSize',text_size,'Color',[0.47 0.67 0.19]);
set(pl1,'Marker','.','LineStyle','none','MarkerEdgeColor','r');
set(pl2,'Marker','.','LineStyle','none','MarkerEdgeColor',[0.47 0.67 0.19]);
set(ax(1),'FontSize',text_size,'YColor','r','Units','pixels', 'Position', [87.5667 67.3000 429.2833 349.9500],'YTick',[0,2.5,5,7.5,10]);
set(ax(2),'FontSize',text_size,'YColor',[0.47 0.67 0.19],'Units','pixels', 'Position', [87.5667 67.3000 429.2833 349.9500],'YTick',[0,2.5,5,7.5,10]*2);
%axis tight;
%print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/bend_inst_11','-dpng','-r0');
%% Bend instability
V = w*56/scale^2*20*1e-18;
conc = zeros(922,1);
num = conc;
for i = 1:922
    df = D(:,i); 
    conc(i) = sum(df(df>0))/(length(df(df>0))*dV);
    num(i) = sum(df);
    %histogram((reshape(atan(results{4,i}./results{3,i}),[],1)),100,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
%    hold on;
end
%plot((1:922)/5,conc);
%plot((1:922)/5,num/V/1e14,'.');
figure('Position',[680 530 600 420]);
[ax,pl1,pl2] = plotyy((1:922)/5,num/V/1e14,VarName1/5,VarName2/scale);
xlabel('Time (s)','FontSize',text_size);
ylabel(ax(1),'{\itc} (\times10^{14} m^{-3})','FontSize',text_size,'Color','b');
ylabel(ax(2),['{\itp}',' (',char(181),'m)'],'FontSize',text_size,'Color','r');
set(pl1,'Marker','.','LineStyle','none','MarkerEdgeColor','b');
set(pl2,'Marker','o','Color','r','LineStyle','-','MarkerEdgeColor','r','MarkerFaceColor','r');
set(ax(1),'FontSize',text_size,'YColor','b','Units','pixels', 'Position', [87.5667 67.3000 429.2833 349.9500]);
set(ax(2),'FontSize',text_size,'YColor','r','Units','pixels', 'Position', [87.5667 67.3000 429.2833 349.9500]);

%set(gca,'FontSize',18); xlim([0 922/5]);
%xlabel('Time (s)','FontSize',18);
%ylabel('{\itc} (\times10^{14} m^{-3})','FontSize',18); box('on');

%ax = gca;
%ax.XTick = [0 pi/2];
%ax.XTickLabel = {'\pi/2','\pi'};859*9
%% Bend instb velocity x distrib
a=1; b = 0; gim = 1; scale = 8192/262.02/10;
quiver(results{1,601}(a:gim:end-b,a:gim:end-b)/scale,...
    (results{2,601}(a:gim:end-b,a:gim:end-b))/scale,...
    results{3,601}(a:gim:end-b,a:gim:end-b),...
    (results{4,601}(a:gim:end-b,a:gim:end-b)),2,'Color','k');
magn=sqrt(results{3,921}.^2+results{4,921}.^2);
%%
per = zeros(921,1);
for i = 921:-1:1
    temp = nanmean(results{4,i});
    %plot(results{1,1}(1,end-length(temp(~isnan(temp)))+1:end),temp(~isnan(temp))); ylim([-10 10]);
    %title(num2str(i));
    %pause(0.1);
    y = temp(~isnan(temp))';
    x = results{1,1}(1,end-length(y)+1:end)';
    
    dx = [1300+i*(200/(921-270)),1800+i*(200/(921-270))];
    
    f = fit(x(x>600),y(x>600),'fourier1');
    per(i) = f.w;
    %plot(x,y); hold on; plot(f);
    %hold off; title(['#=',num2str(i),' period = ',num2str(2*pi/per(i))]);
    %pause(0.1);
    
    %f = fit(results{1,1}(1,end-100:end)',temp(1,end-100:end)','fourier1');
    %plot(results{1,1}(1,end-100:end)',temp(1,end-100:end)'); hold on;
end

%% 3D plot Bend instability
load('/media/turbik/DATA/PhD/DATA/2017_Bacteria_c,s,v-stripes/NewRM_cells/Cell_fluoresc_5fps/PTV_all_7sigm_0.5_16area.mat');
%%
v = zeros(921,2);
for i = 1:921
    temp = sqrt(nanmean(results{3,i}).^2 + nanmean(results{4,i}).^2);
    hold on;
    %plot3((i+ones(length(temp(~isnan(temp))),1))/5,700-(results{1,1}(1,end-length(temp(~isnan(temp)))+1:end))/scale,temp(~isnan(temp)),'LineWidth',1); 
    %ylim([-10 10]);
    %title(num2str(i));
    v(i,1) = nanmean(nanmean(results{3,i}));
    v(i,2) = nanmean(nanmean(results{4,i}));
end
grid on; view(-85,45);
xlabel(gca, 'Time (s)'); ylabel(['{\itx} ',' (',char(181),'m)']); zlabel(['{\itv_y}',' (',char(181),'m/s)']);
set(gca,'FontSize',text_size);

%%
plot3((i+ones(length(temp(~isnan(temp))),1))/5,700-(results{1,1}(1,end-length(temp(~isnan(temp)))+1:end))/scale,temp(~isnan(temp)),'LineWidth',1);
%% 0.082853*(1:921)'+85.693
dV = 2100*20*160/scale*10^(-18);

%x = (0.082853*(1:921)'+85.693)/dV; y = smooth(2*pi./smooth(per)/scale,10);

%plot(x(240:end), smooth(y(240:end),1),'.-')
plot(dat(239:end,1)/dV,smooth(dat(239:end,2),1)/1e6,'.'); hold on;
plot(dat(239:500,1)/dV,(140000./(dat(239:500,1)/dV - 4.875e13)).^0.5,'r-');
text_size=18;

xlabel('Bacteria count','FontSize',text_size);
ylabel(['{\it\Lambda} ',' (',char(181),'m)'],'FontSize',text_size);
set(gca,'FontSize',text_size); ylim([0 0.001]);
%print('/media/turbik/DATA/PhD/PUBLICATIONs/17_Bacteria_c-stripes/per_vs_bact','-dpng','-r0');
%axis tight;
%%
load('d.mat');
text_size=18;
figure('Position',[680 530 640 440]);
[ax,pl1,pl2] = plotyy(dat(239:end,1)/dV,smooth(dat(239:end,2),1)/1e6,dat(239:end,1)/dV,smooth(test(239:end-1,2)/scale/0.2,5));
%0.082853*(1:921)'+85.693,smooth(2*pi./smooth(per)/scale,10),0.082853*(1:921)'+85.693,smooth(test(1:end-1,2)/scale/0.2,5));
%axis tight;
hold(ax(1),'on');
plot(dat(239:500,1)/dV,(140000./(dat(239:500,1)/dV - 4.875e13)).^0.5,'r-','LineWidth',1);
%plot(0.082853*(-100:921)'+85.693,((3.06549e-6./((0.082853*(-100:921)'+85.693)- 79.1748)).^0.5)*1e6,'k--');

xlabel('{\itc} (\times 10^{13} m^{-3})','FontSize',text_size);
ylabel(ax(1),['{\it\lambda} ',' (',char(181),'m)'],'FontSize',text_size);
ylabel(ax(2),['\langle{\itv}\rangle',' (',char(181),'m/s)'],'FontSize',text_size,'Color','r');
set(pl1,'Marker','.','LineStyle','none','MarkerEdgeColor','b');
set(pl2,'Marker','.','LineStyle','-','MarkerEdgeColor','m');
set(ax(1),'FontSize',text_size,'YColor','b','Units','pixels', 'Position', [117.5667 80 429.2833 349.9500],'YLim',[0 1.2e-3],'YTick',(0:0.2:1.2)*1e-3,'XTick',(5:8)*1e13,'XTickLabel',5:8,'YTickLabel',0:200:1200);
set(ax(2),'FontSize',text_size,'YColor','m','Units','pixels', 'Position', [117.5667 80 429.2833 349.9500],'YLim',[0,15],'YTick',[0,2.5,5,7.5,10]*2);
%% ^2 of concentration
figure('Position',[680 530 640 420]);
[ax,pl1,pl2] = plotyy(((0.082853*(1:921)'+85.693)/2.1484e-13/6e13).^2,smooth(2*pi./smooth(per)/scale,10),0.082853*(1:921)'+85.693,smooth(test(1:end-1,2)/scale/0.2,5));
hold on;
plot(3.06549e-6/(x - 79.1748))^0.5);
%axis tight;
xlabel('Bacteria count','FontSize',text_size);
ylabel(ax(1),['{\it\lambda} ',' (',char(181),'m)'],'FontSize',text_size);
ylabel(ax(2),['\langle{\itv}\rangle',' (',char(181),'m/s)'],'FontSize',text_size,'Color','r');
set(pl1,'Marker','.','LineStyle','-','MarkerEdgeColor','b');
set(pl2,'Marker','.','LineStyle','-','MarkerEdgeColor','r');
set(ax(1),'FontSize',text_size,'YColor','b','Units','pixels', 'Position', [117.5667 67.3000 429.2833 349.9500],'YLim',[0,1100]);
set(ax(2),'FontSize',text_size,'YColor','r','Units','pixels', 'Position', [117.5667 67.3000 429.2833 349.9500],'YLim',[0,20],'YTick',[0,2.5,5,7.5,10]*2);

%% PolScope
P = cell(5,2);
tStep = [0.2,0.5,1,2,5];
for i=1:5
    if i<=2, fps = 6;
    else fps = 3; end
    if i<=4, key = 1;
    else key = 1; end
    P{i,1} = tStep(i);
    P{i,2}(:,2) = temp(:,1);
    P{i,2}(:,1) = 36-tStep(i)/60*fps*(0:length(temp)-1);
end
% P
figure;
for i=1:5
    subplot(2,1,1);
    plot(P{i,2}(:,1),(P{i,2}(:,2))); hold on;
end
legend('0.2','0.5','1','2','5');xlabel('Temperature (^oC)'); ylabel('Optical phase retardation (nm)'); 
title('PolScope');
xlim([30 37]);

%% S
for i=1:5
    subplot(2,1,2);
    plot(S{i,2}(:,1),smooth(S{i,2}(:,2))); hold on;
end
legend('0.2','0.5','1','2','5');
xlim([30 37]); ylim([0 1]);
xlabel('Temperature (^oC)'); ylabel('Order parameter'); title('POM');
%plot(S);
%% Particle Transport by bacteria
d1 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/48.24s.csv',',',1,1);
d2 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/2.csv',',',1,1);
d3 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/bent_inst_1.csv',',',1,1);
sc = 8192/185.2/10;

figure('Position',[600 500 600 400],'PaperPositionMode','auto','Color','w');
%set(gca,'FontSize',text_size,'YLim',[-5 12],'XTick',0:250:1000,...
%    'Unit','pixel','Position',[125 100 450 250],'XMinorTick','on','YMinorTick','on'); 
plot(d1(:,3)/sc,d1(:,4)/sc-91,'bo',d2(:,3)/sc,d2(:,4)/sc-91,'r*','LineWidth',1.5); axis equal;% ,d3(:,3)/sc,d3(:,4)/sc,'^'

data5g_1 = [d1(:,3)/sc,d1(:,4)/sc-91];
data5g_2 = [d2(:,3)/sc,d2(:,4)/sc-91];

data5g_1 = data5g_1(end:-1:1,:);
data5g_2 = data5g_2(end:-1:1,:);

set(gca,'FontSize',18,'Unit','pixel','Position',[125 100 450 300],'LineWidth',1.5);%,'Position',[0.1300 0.2227 0.7750 0.
ylabel(['{\ity} (',char(181),'m)']);
xlabel(['{\itx} (',char(181),'m)']);
l = legend({'Single colloid','Multiple colloids'},'Box','off');% ,'Single colloid, undulating jet'
l.FontSize = 14;
pause(2);
print('17_Bacteria_c-stripes/trasnport_coordinates','-dpng','-r300'); close(gcf);
%% data into csv
fname = 'file_data_5g.csv';
fid = fopen(fname,'w'); 
fprintf(fid,'Single colloid:\nx(um),y(um)\n'); fclose(fid);
dlmwrite(fname,data5g_1,'-append');

fid = fopen(fname,'a'); 
fprintf(fid,'\nMultiple colloids:\nx(um),y(um)\n'); fclose(fid);
dlmwrite(fname,data5g_2,'-append');
%% transport velocity
figure('Position',[600 500 600 400],'PaperPositionMode','auto','Color','w');
tstep = 10/8.5;
plot((1:length(d1)-1)/tstep,-smooth(diff(sqrt(d1(:,4).^2+d1(:,3).^2)/sc/tstep)),'bo-',...
    (1:length(d2)-1)/tstep,-smooth(diff(sqrt(d2(:,4).^2+d2(:,3).^2)/sc/tstep)),'r*-','LineWidth',1.5); 
%,...    (1:length(d3)-1)/tstep,-smooth(diff(sqrt(d3(:,4).^2+d3(:,3).^2)/sc/tstep)),'^-'    

data5h_1 = [(1:length(d1)-1)'/tstep,(-smooth(diff(sqrt(d1(:,4).^2+d1(:,3).^2)/sc/tstep)))];
data5h_2 = [(1:length(d2)-1)'/tstep,(-smooth(diff(sqrt(d2(:,4).^2+d2(:,3).^2)/sc/tstep)))];

%xlabel('t (s)'); ylabel(['v ',' (',char(181),'m/s)']);
set(gca,'Unit','pixel','Position',[125 100 450 300],'FontSize',18,'LineWidth',1.5);%,'Position',[0.1300 0.2227 0.7750 0.
ylabel(['{\itv}_{c,x} (',char(181),'m/s)']);
xlabel('{\itt} (s)');

l = legend({'Single colloid','Multiple colloids'},'Box','off');
l.FontSize = 14;
ylim([0 12]);
mvd1 = mean(diff(sqrt(d1(:,4).^2+d1(:,3).^2)/sc/tstep));
mvd2 = mean(diff(sqrt(d2(:,4).^2+d2(:,3).^2)/sc/tstep));
mvd3 = mean(diff(sqrt(d3(:,4).^2+d3(:,3).^2)/sc/tstep));
pause(2);
print('17_Bacteria_c-stripes/trasnport_velocity','-dpng','-r0'); close(gcf);

%% data into csv
fname = 'file_data_5h.csv';
fid = fopen(fname,'w'); 
fprintf(fid,'Single colloid:\nt(s),v_{c,x}(um/s)\n'); fclose(fid);
dlmwrite(fname,data5h_1,'-append');

fid = fopen(fname,'a'); 
fprintf(fid,'\nMultiple colloids:\nt(s),v_{c,x}(um/s)\n'); fclose(fid);
dlmwrite(fname,data5h_2,'-append');
%% Diffusion / Particle Transport by bacteria
d1 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/48.24s.csv',',',1,1);
d2 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/2.csv',',',1,1);
d3 = dlmread('2017_Bacteria_c,s,v-stripes/Particle_transport/bent_inst_1.csv',',',1,1);
sc = 8192/185.2/10;
%set(gca,'FontSize',text_size,'YLim',[-5 12],'XTick',0:250:1000,...
%    'Unit','pixel','Position',[125 100 450 250],'XMinorTick','on','YMinorTick','on'); 

x = d1(:,3)/sc;
y = d1(:,4)/sc-91;
[dX2,dY2,DT] = MSD_funct(x,y,10/8.5);
f1 = fit(DT,dX2*1e-12,'power1','StartPoint',[100e-12, 1.581]);
f2 = fit(DT,dX2*1e-12,'a*x^2+b*x');
f3 = fit(DT,dX2*1e-12,'a*(x-x*exp(-x/b))');

%%
figure('Position',[600 500 600 400],'PaperPositionMode','auto','Color','w');
p2 = plot((0.1:0.1:31),f1(0.1:0.1:31)); hold on;
p1 = plot(DT,dX2*1e-12,'bo',DT,dY2*1e-12,'rsq','LineWidth',1.5); %hold on; 
data5i = [DT,dX2*1e-12,dY2*1e-12];
set(p1(1),'MarkerFaceColor','b');
set(p1(2),'MarkerFaceColor','r');
xlim([0 35]);
ylim([0 3e-8]);
set(p2,'Color','k','LineWidth',1.5,'LineStyle','-');
set(gca,'FontSize',18,'Unit','pixel','Position',[125 100 450 300],'LineWidth',1.5);%,'Position',[0.1300 0.2227 0.7750 0.
set(gca,'YTickLabel',get(gca,'YTick')/1e-8);
ylabel(['MSD (\times10^{-4} ',char(181),'m^2/s)']);
xlabel('\Delta{\itt} (s)');
l = legend([p1(1),p1(2),p2],{'\langle\Delta{\itx}^2\rangle','\langle\Delta{\ity}^2\rangle','Fitting'},'Box','off','Location','southeast');% ,'Single colloid, undulating jet'
l.FontSize = 14;

%pause(2);
h = axes('Unit','pixel','Position',[205 250 166 135]);
loglog(h,DT,dX2*1e-12,'bo','LineWidth',1.5,'MarkerFaceColor','b'); hold on;
loglog(h,DT,dY2*1e-12,'rsq','LineWidth',1.5,'MarkerFaceColor','r');
%set(p1(1),'MarkerFaceColor','b');
ylim(h,[1e-13 1e-7]);
%ylabel(h,['MSD (',char(181),'m^2/s)']);
ylabel(h,['MSD (m^2/s)']);
xlabel(h,'\Delta{\itt} (s)');
set(h,'LineWidth',1.5,'FontSize',12,'TickLength',[0.025 0],'YTick',[1e-12,1e-10,1e-8]);
%loglog(h,[1e1 10^(1.5)],[1e-10 10^(-9)],'k','LineWidth',2); hold on;
%loglog(h,[10^1 10^1.5],[1e-10 1e-10],'k','LineWidth',2);
annotation('textbox',[0.35 0.68,.2 .2],'String','\propto\Delta{\itt}^2','FitBoxToText','on','LineStyle','none');
annotation('textbox',[0.5 0.62,.2 .2],'String','\propto\Delta{\itt}^{1.5}','FitBoxToText','on','LineStyle','none');
loglog(h,[10^0.5 10^(1.5)],[1e-10 10^(-8.5)],'k','LineWidth',2); hold on;
loglog(h,[10^0 10^(1)],[10^(-9.5) 10^(-7.5)],'k','LineWidth',2);
%loglog(h,[10^1 10^1.5],[1e-10 1e-10],'k','LineWidth',2);
pause(2);
print('17_Bacteria_c-stripes/trasnport_coordinates_3','-dpng','-r300'); close(gcf);
%% data for csv
fname = 'file_data_5i.csv';
fid = fopen(fname,'w'); 
fprintf(fid,'Delta t(s),MSD_x(m^2/s),MSD_y(m^2/s)\n'); 
fclose(fid);
dlmwrite(fname,data5i,'-append');

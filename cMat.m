function []=cMat(y_test,y_predict)
%% confusion matrix
 Tu=find(y_test==0);
 Pu=find(y_predict==0);
 Td=find(y_test==1);
 Pd=find(y_predict==1);

 Ltp=length(find(ismember(Tu,Pu)));
 Lfn=length(find(ismember(Tu,Pd)));
 Lfp=length(find(ismember(Td,Pu)));
 Ltn=length(find(ismember(Td,Pd)));

 CM=[Ltp Lfn;Lfp Ltn];
ac_cm=(Ltp+Ltn)/(Ltp+Lfn+Lfp+Ltn)
%% Plot
u=length(Tu);
d=length(Td);
CM_n=[CM(1,1)/u CM(1,2)/u;CM(2,1)/d, CM(2,2)/d] %normalized confusion matrix
labels = {'Undamaged', 'Damaged'};
figure
 % plotting the colors
imagesc(CM_n);
%title('Confusion Matrix');
ylabel('True Class'); xlabel('Predicted Class');
set(gca,'FontSize',16)
% set the colormap
load('colormapsavefile.mat');
colormap(myColormap);
colorbar
% Create strings from the matrix values and remove spaces
textStrings = num2str([CM_n(:), CM(:)], '%.2f%\n%d\n');
textStrings = strtrim(cellstr(textStrings));

% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:2);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center', 'Fontsize',14);

% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));

% Choose white or black for the text color of the strings so
% % they can be easily seen over the background color
% textColors = repmat(CM_n(:) > midValue,1,3);
% set(hStrings,{'Color'},num2cell(textColors,2));
numlabels=2;
% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'TickLength',[0 0]);
end

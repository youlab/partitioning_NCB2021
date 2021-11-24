% paper color scheme

function pColor = paperColor
pColor = [156, 133, 192; 113, 96, 140; 
    165,181,146; 140, 150, 127;
    231, 188, 41; 170, 137, 27; 
    243, 164, 71; 179, 119, 50;
    128, 158, 194; 95, 117, 145;
    ]/255;

% visualize the colors
figure(345)

P=bar(ones(2,10),'stacked');
for i = 1:10
    set(P(i),'FaceColor',pColor(i,:));
end
filename = './ZDT1/ZDT1.txt';
[x,y]=textread(filename,'%f %f');

scatter(x,y,'k');

hold on;

filename = './ZDT1/my.txt';
[x,y]=textread(filename,'%f %f');

scatter(x,y,'g');
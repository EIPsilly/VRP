filename = './ZDT4/ZDT4.txt';
[x,y]=textread(filename,'%f %f');

scatter(x,y,'k');

hold on;

filename = './ZDT4/my.txt';
[x,y]=textread(filename,'%f %f');

scatter(x,y,'r');
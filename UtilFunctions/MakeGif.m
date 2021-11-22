function [outputArg1,outputArg2] = MakeGif(fname,dirTF)
if nargin<1
    fname = 'youforgettogivemeaname.gif';
end
if nargin<2
    dirTF = '.\';
end

toImg = dir([dirTF,'*).png']);

names = {};
for k = 1:length(toImg)
    names{k} = toImg(k).name;
end
names = [names{1},names]

for n = 1:length(names)
    I = imread(names{n});
    
    [imind,cm] = rgb2ind(I,256);
    
    if n == 1
        imwrite(imind,cm,fname,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',1);
    end
end





end


function [en, ecropped, eseamed] = gradient_sal_local()
close all;
%-------Read Image File
imname = input('Enter the name of the image: ');
I = imread(imname);
img = I;
original_image = I;
Im = rgb2gray(I);
s = size(I);
disp(s);
%-----------------Y contains the destination location for output images.
Y = 'C:\Users\Arushi\Desktop\Files\DIP\';
%d = input('Enter 0 for vertical, 1 for horizontal,2 for vertical followed by horizontal, 3 for horizontal followed by vertical: ');
%---------right now only for vertical seam carving.
d = 0;
if ((d == 0) || (d == 1))
    %inp = input('No of rows/columns you want to remove: ');
    inp = input('No of columns you want to remove: ');
else
    inp1 = input('No of columns you want to remove: ');
    inp2 = input('No of rows you want to remove: ');
end

%----------Input Parameters
small = input('small block size: ');
big = input('big block size: ');
ef = input('edge factor: ');
df = input('ID factor: ');

%----------Calculates gradient map
        Im = rgb2gray(I);
        eI = find_gradient(Im);
        imshow(eI,[]);
        X = strcat('gradientmap', imname);
        file = fullfile(Y,X);
        saveas(gcf, file);

%-----------Calculates global saliency        
        ecolor = global_color(I);
        figure(2);imshow(ecolor,[]);
        X = strcat('global_sal_', imname);
        file = fullfile(Y,X);
        saveas(gcf, file);

%----------Calculates local saliency for 2 block sizes 
        localS = zeros(size(ecolor));
        localS2 = localS;
        if small
        localS = local_sal(I, small);
        localS = localS./max(max(localS))*max(max(ecolor));
        figure(3);imshow(localS,[]);
        X = strcat('localsmall_sal_', imname);
        file = fullfile(Y,X);
        saveas(gcf, file);
        end;
        if big       
        localS2 = local_sal(I, big);
        localS2 = localS2./max(max(localS2))*max(max(ecolor));
        figure(4);imshow(localS2,[]);
        X = strcat('localbig_sal_', imname);
        file = fullfile(Y,X);
        saveas(gcf, file);
        end
  
%-------------adds gradient, global and local saliency for final energy map
        ecolor = ecolor + (localS2 + localS)*2;
        figure(5);imshow(ecolor,[]);
        X = strcat('added_sal_', imname);
        file = fullfile(Y,X);
        saveas(gcf, file);
        
sz = size(eI);
eI = eI./max(max(eI))*max(max(ecolor));
en = eI(2:sz(1)-1, 2:sz(2)-1)*ef + ecolor;
figure(6);imshow(en,[]);
X = strcat('total_energy_', imname);
file = fullfile(Y,X);
saveas(gcf, file);

%----------Calulating optimal seam using Dynamic Programming
%-------MV = vertically summped up energy map
MV = vertical_seam_energy(en,s);
ecropped = en(:,inp/2+1 : sz(2) - inp/2);
eseamed = en;

if (d == 0 || d == 1)
for loop=1:inp
     s = size(I);
     if (d == 0)
        
%         MV = vertical_seam_energy(en,s);
        seam = locate_optimal_vertical_seam(MV,s);
        [J, NV] = remove_optimal_vertical_seam(I,MV, seam, s,df);
        [eseamed] = remove_optimal_vertical_energy_seam(eseamed, seam, size(eseamed));
        MV = NV;
%         loop = loop+1;
     else
        e1 = find_gradient(I(:,:,1));
        e2 = find_gradient(I(:,:,2));
        e3 = find_gradient(I(:,:,3));
        eI = e1+e2+e3;
        en = find_entropy(eI,rgb2gray(I));
        eI = (en);
        MH = horizontal_seam_energy(eI,s);
        seam = locate_optimal_horizontal_seam(MH,s);
        [J, NH] = remove_optimal_horizontal_seam(I,MH, seam, s);        
%         MH = NH;
     end
     I = J;
end
end

disp(size(original_image));
S = size(I);
disp(S);

figure()
imshow(original_image);
figure()
imshow(I);
X = strcat('final_', imname);
file = fullfile(Y,X);
imwrite(I,file);
end

%................Gradient Filters...................%
function [eI] = find_gradient(Im)
Im = imfilter(Im, fspecial('gaussian', 5, 5), 'symmetric', 'conv');

gx = [-1,0,1; -2,0,2; -1,0,1];
gy = [-1,-2,-1; 0,0,0; 1,2,1];

ex = abs(conv2(gx,double(Im)));
ey = abs(conv2(gy,double(Im)));

eI = ex+ey;
mn = mean(mean(eI));
eI(eI < mn) = 0;
end

function [sm] = global_color(img)
gfrgb = imfilter(img, fspecial('gaussian', 3, 3), 'symmetric', 'conv');

cform = makecform('srgb2lab');
lab = applycform(gfrgb,cform);

l = double(lab(:,:,1)); lm = mean(mean(l));
a = double(lab(:,:,2)); am = mean(mean(a));
b = double(lab(:,:,3)); bm = mean(mean(b));

sm = (l-lm).^2 + (a-am).^2 + (b-bm).^2;
% figure();imshow(sm,[]);
end

function [temp] = local_sal(Im, inp)
% Im = imread('d.jpg');
[s1, s2, layers] = size(Im);
% inp = 20;
temp = zeros(s1,s2);

for x = 1:inp:s1
    if(x+inp-1 < s1)
        for y = 1:inp:s2
%             disp(s1);
%             disp(s2);
%             disp(x);
%             disp(y);
%             disp('************************');
            if(y+inp-1 < s2)                           
                temp(x:x+inp-1, y:y+inp-1) = global_color(Im(x:x+inp-1, y:y+inp-1,:));
%                 en(x:x+inp-1, y:y+inp-1)= en(x:x+8, y:y+inp-1) +temp;
            end
        end
        if y<s2
            temp(x:x+inp-1, y:s2) = global_color(Im(x:x+inp-1, y:s2,:));
%             en(x:x+inp-1, y:s2)= en(x:x+inp-1, y:s2) +temp;
        end
    end
end
if x < s1
    for y = 1:inp:s2
        if(y+inp-1 < s2)
        temp(x:s1, y:y+inp-1) = global_color(Im(x:s1, y:y+inp-1,:));
%         en(x:s1, y:y+inp-1)= en(x:s1, y:y+inp-1) +temp;
        end
    end
    if y<s2
        temp(x:s1, y:s2) = global_color(Im(x:s1, y:s2,:));
%         en(x:s1, y:s2)= en(x:s1, y:s2) +temp;
    end
end

end

function [en] = find_entropy(eI, Im)
[s1, s2] = size(Im);
en = eI;
for x = 1:9:s1
    if(x+8 < s1)
        for y = 1:9:s2
            if(y+8 < s2)
                temp = entropy(Im(x:x+8, y:y+8));
                en(x:x+8, y:y+8)= en(x:x+8, y:y+8) +temp;
            end
        end
        if y<s2
            temp= entropy(Im(x:x+8, y:s2));
            en(x:x+8, y:s2)= en(x:x+8, y:s2) +temp;
        end
    end
end
if x < s1
    for y = 1:9:s2
        if(y+8 < s2)
        temp = entropy(Im(x:s1, y:y+8));
        en(x:s1, y:y+8)= en(x:s1, y:y+8) +temp;
        end
    end
    if y<s2
        temp= entropy(Im(x:s1, y:s2));
        en(x:s1, y:s2)= en(x:s1, y:s2) +temp;
    end
end
end

%.............Calculating Vertical Seam Energies ..............%
function [MV] = vertical_seam_energy(eI, s)
MV(1,:) = eI(1,:);
for i=2:s(1)
    for j = 1:s(2)
        if(j == 1)
            MV(i,j)  = eI(i,j) + min(MV(i-1,j),MV(i-1,j+1));
        elseif(j == s(2))
            MV(i,j)  = eI(i,j) + min(MV(i-1,j-1),MV(i-1,j));
        else
            MV(i,j)  = eI(i,j) + min(MV(i-1,j-1),min(MV(i-1,j),MV(i-1,j+1)));
        end
    end
end


% m = min(min(MV));
% x = max(max(MV));

MV = MV(:,1:(s(2)-1));
% imshow(MV, [m x]);

end

function [MH] = horizontal_seam_energy(eI, s)
MH(:,1) = eI(:,1);
for i=1:s(1)
    for j = 2:s(2)
        if(i == 1)
            MH(i,j)  = eI(i,j) + min(MH(i,j-1),MH(i+1,j-1));
        elseif(i == s(1))
            MH(i,j)  = eI(i,j) + min(MH(i-1,j-1),MH(i,j-1));
        else
            MH(i,j)  = eI(i,j) + min(MH(i-1,j-1),min(MH(i,j-1),MH(i+1,j-1)));
        end
    end
end


% m = min(min(MH));
% x = max(max(MH));

MH = MH(1:(s(1)-1),:);
% imshow(MH, [m x]);

end

function [seam] = locate_optimal_vertical_seam(MV,s)
    clear J;
    k = min(MV(s(1),:));
    f = find(MV(s(1),:) == k);
    f = f(1);
    seam = zeros(s(1),1);
    z = s(1);
    seam(z) = f;
    while(z ~= 1)
        z = z - 1;
        small = MV(z,f);
        if (f == 1)
            if (MV(z,f+1) < small)
                f = f + 1;
            end
        elseif (f == (s(2)-1))
            if(MV(z,f-1) < small)
                f = f - 1;
            end
        else
            if(MV(z,f-1) < small)
                f = f - 1;
            elseif (MV(z,f+1) < small)
                f = f + 1;
            end
        end
        seam(z) = f;
    end
end

% function [seam,num] = locate_optimal_vertical_seam(MV,s)
%     clear J;
%     k = min(MV(s(1),:));
%     F = find(MV(s(1),:) == k);
%     num = size(F);   
%     seam = zeros(s(1),num);
%     z = s(1);
%     for x= 1:num
%     f = F(x);
%     seam(z,x) = f;
%     while(z ~= 1)
%         z = z - 1;
%         small = MV(z,f);
%         if (f == 1)
%             if (MV(z,f+1) < small)
%                 f = f + 1;
%             end
%         elseif (f == (s(2)-1))
%             if(MV(z,f-1) < small)
%                 f = f - 1;
%             end
%         else
%             if(MV(z,f-1) < small)
%                 f = f - 1;
%             elseif (MV(z,f+1) < small)
%                 f = f + 1;
%             end
%         end
%         seam(z,x) = f;
%     end
%     end
% end

function [seam] = locate_optimal_horizontal_seam(MH,s)
    clear J;
    k = min(MH(:,s(2)));
    f = find(MH(:,s(2)) == k);
    f = f(1);
    seam = zeros(1,s(2));
    z = s(2);
    seam(1,z) = f;
    while(z ~= 1)
        z = z - 1;
        small = MH(f,z);
        if (f == s(1)-1)
            if(MH(f-1,z) < small)
                f = f - 1;
            end
        elseif (f == 1)
            if(MH(f+1,z) < small)
                f = f + 1;
            end
        else
            if(MH(f-1,z) < small)
                f = f - 1;
            elseif (MH(f+1,z) < small)
                f = f + 1;
            end
        end
        
        seam(1,z) = f;
    end
end

function [J, NV] = remove_optimal_vertical_seam(I, MV, seam, s,df)
    for i=1:s(1)
        J(i,1:seam(i)-1,:) = I(i,1:seam(i)-1,:);
        NV(i,1:seam(i)-2) = MV(i,1:seam(i)-2);
        if(seam(i) - 1)
        NV(i,seam(i)-1) = MV(i,seam(i)-1)+ df*MV(i, seam(i));
        end
        J(i,seam(i):(s(2)-1),:) = I(i,(seam(i)+1):s(2),:);
        if(seam(i)+1 < s(2))
        NV(i,seam(i)) = MV(i,seam(i)+1)+ df*MV(i, seam(i));
        end
        NV(i,seam(i)+1:(s(2)-2)) = MV(i,(seam(i)+2):s(2)-1);
    end
end

function [J] = remove_optimal_vertical_energy_seam(I, seam, s)
    for i=1:s(1)
        J(i,1:seam(i)-1) = I(i,1:seam(i)-1);        
        J(i,seam(i):(s(2)-1),:) = I(i,(seam(i)+1):s(2),:);        
    end
end


function [J, NH] = remove_optimal_horizontal_seam(I, MH, seam, s)
    for i=1:s(2)
        J(1:seam(i)-1,i,:) = I(1:seam(i)-1,i,:);
        NH(1:seam(i)-1,i) = MH(1:seam(i)-1,i);
        J(seam(i):(s(1)-1),i,:) = I((seam(i)+1):s(1),i,:);
        NH(seam(i):(s(1)-2),i) = MH((seam(i)+1):s(1)-1,i);
    end
end

function [] = seam_carving_insertion()

%--------------Read Image and parameters
imname = input('Enter the name of the image: ');
%lessname = input('Enter the name of the less image: ');
%seamname = input('Enter the name of the image seams: ');
% I = imrotate(imread(imname), -90);
I = (imread(imname));
img = I;
original_image = I;
Im = rgb2gray(I);
bigI = I;
s = size(I);
disp(s);
d = input('Enter 0 for vertical, 1 for horizontal seams: ');

    inp = input('No of rows/columns you want to remove: ');

%---------finding image gradients and entropy    
eI2 = find_gradient(Im);
en2 = find_entropy(eI2, Im);
% en = find_entropy(eI, Im);
% eI = (en);
% e2 = find_gradient(I(:,:,2));
% e3 = find_gradient(I(:,:,3));
% eI = e1+e2+e3;
% if(d == 0)
%      MV = vertical_seam_energy(eI,s);
% elseif(d == 1)
%     MH = horizontal_seam_energy(eI,s);
% else
%     MV = vertical_seam_energy(eI,s);
%     MH = horizontal_seam_energy(eI,s);
% end

for loop=1:inp
     s = size(I);
     if (d == 0)
         %------------image gradient and energy map calculate after every
         %seam removal - good results but time consuming.
        e1 = find_gradient(I(:,:,1));
        e2 = find_gradient(I(:,:,2));
        e3 = find_gradient(I(:,:,3));
        eI = e1+e2+e3;
%         en = HoG(eI);
%         eI = (en);
        MV = vertical_seam_energy(eI,s);
        seam(:, loop) = locate_optimal_vertical_seam(MV,s);
        [J, NV] = remove_optimal_vertical_seam(I,MV, seam(:,loop), s);
%         MV = NV;
%         loop = loop+1;
     else
         %------------image gradient and energy map calculate after every
         %seam removal - good results but time consuming.
        e1 = find_gradient(I(:,:,1));
        e2 = find_gradient(I(:,:,2));
        e3 = find_gradient(I(:,:,3));
        eI = e1+e2+e3;
        en = find_entropy(eI,rgb2gray(I));
        eI = (en);
        MH = horizontal_seam_energy(eI,s);
        seam(loop, :) = locate_optimal_horizontal_seam(MH,s);
        [J, NH] = remove_optimal_horizontal_seam(I,MH, seam(loop, :), s);
%         MH = NH;
     end
     I = J;
end

[bigI, newI] = insert_vertical_seam(bigI, seam, inp);

disp(size(original_image));
S = size(I);
disp(S);
%------------displaying results
figure();
imshow(original_image);
title('original image');
figure();
imshow(I);
title('seam_carved image');
% imwrite(imrotate(I, 90),lessname); 
%imwrite(I,lessname); 
figure();
imshow(newI);
title('seams');
% imwrite(imrotate(newI, 90),seamname); 
%imwrite(newI,seamname);
figure();
imshow(bigI);
title('enlarged image');
end

%................Gradient Filters...................%
function [eI] = find_gradient(Im)
gx = [-1,0,1; -2,0,2; -1,0,1];
gy = [-1,-2,-1; 0,0,0; 1,2,1];

ex = abs(conv2(gx,double(Im)));
ey = abs(conv2(gy,double(Im)));

eI = ex+ey;
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

function [eHoG] = HoG(eI)
[s1, s2] = size(eI);
% eHoG = eI;
for x = 1:11:s1
    if(x+10 < s1)
        for y = 1:11:s2
            if(y+10 < s2)
                temp = hist(eI(x:x+10, y:y+10),8);
                temp = max(sum(temp,2));
                eHoG(x:x+10, y:y+10)= eI(x:x+10, y:y+10)./temp;
            end
        end
        if y<s2
            temp= hist(eI(x:x+10, y:s2),8);
            temp = max(sum(temp,2));
            eHoG(x:x+10, y:s2)= eI(x:x+10, y:s2)./temp;
        end
    end
end
if x < s1
    for y = 1:11:s2
        if(y+10 < s2)
        temp = hist(eI(x:s1, y:y+10),8);
        temp = max(sum(temp,2));
        eHoG(x:s1, y:y+10)= eI(x:s1, y:y+10)./temp;
        end
    end
    if y<s2
        temp= hist(eI(x:s1, y:s2),8);
        temp = max(sum(temp,2));
        eHoG(x:s1, y:s2)= eI(x:s1, y:s2)./temp;
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

function [J, NV] = remove_optimal_vertical_seam(I, MV, seam, s)
    for i=1:s(1)
        J(i,1:seam(i)-1,:) = I(i,1:seam(i)-1,:);
        NV(i,1:seam(i)-1) = MV(i,1:seam(i)-1);
        J(i,seam(i):(s(2)-1),:) = I(i,(seam(i)+1):s(2),:);
        NV(i,seam(i):(s(2)-2)) = MV(i,(seam(i)+1):s(2)-1);
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

function [bigI, newI] = insert_vertical_seam(I, seam, inp)
s = size(I);
newI = I;
temp = zeros(s);
for x= 1:inp
    for y = 1:x-1
       temp(:, x) = temp(:,x) + double(seam(:, x) >= seam(:, x-y));
    end
    seam(:, x)= seam(:, x) + temp(:, x);
end
for x= 1:inp
    for y = 1: s(1)
        newI(y,seam(y,x), :) =0;
        newI(y,seam(y,x), 1) =255;
    end
    s = size(I);
    for i=1:s(1)
        bigI(i,1:seam(i, x),:) = I(i,1:seam(i, x),:);
        bigI(i, seam(i, x)+1,:) = I(i,seam(i, x),:);
        bigI(i,seam(i, x)+2:(s(2)+1),:) = I(i,(seam(i, x)+1):s(2),:);
    end
    I = bigI;
end
end
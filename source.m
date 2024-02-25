fprintf("Enter the coordinate value of center.\n");
a = input('x: ');
b = input('y: ');
r = input('Enter the radius: ');
dt = input('Enter the number of detectors & Sources (in Odd num): ');    %is dt always even?
delta = input("Enter distance between detectors ");
SDD = input("Enter SDD: ");
SCD = input("Enter SCD: ");

view = input("How many views?: ");
IN = input("Enter Attenuation coefficient of circle: ");
dsxArr = [-((dt-1)*delta/2):delta:((dt-1)*delta/2)];  % x values of Detector & Source in Cart esian coordinate
angch = 2*pi/(view);
FPres = zeros(view, dt); % Store Forward Projection Result.

cnt=0;

%Forward Projection
%Step 1. Find projection from each angle. 
for t = 0:angch:2*pi-angch % rotate angle
    A = [cos(t) sin(t) ; -sin(t) cos(t)];
    cnt = cnt+1;
    for i = 1:dt   % #detector
        q1 = A*[dsxArr(i); b - SCD];
        q2 = A*[dsxArr(i); b - SCD + SDD];
        FPres(cnt, i) = IN*FP(a,b,r,q1,q2);
    end
end
%sinogram_size = size(FPres);
subplot(1, 3, 1);   
imshow(FPres, [0,0.005]);
title('Original Sinogram');

pad_size = 2 ^ ceil(log2(2*dt)); 
FPres = fft(FPres, pad_size, 2);

% Step 2. h(n*theta)를 구현 후 Fourier transform
rampFilter = zeros(1, pad_size/2+1);   %rampFilter
for row = 1:(pad_size/2+1)   %1-base indexing이라서 1칸씩 오른쪽으로 shift
    if(row==1)
        rampFilter(row) = 1/(4 * (delta)^2);
    elseif (mod(row,2)==0)
        rampFilter(row) = -1/(((row-1)*pi*delta)^2);
    end
end
% Size 조정: 완성한 Filter을 뒤집어서 반대방향으로도 붙여준다.
reverserampFilter = fliplr(rampFilter);
reverserampFilter(:,pad_size/2+1) = [];
rampFilter = [reverserampFilter rampFilter];
rampFilter = fft(rampFilter, pad_size); %이 fft(X,n)을 하면 푸리에 변환 시, 함수 안에서 자체적으로 X를 n만큼 늘린 값으로 푸리에변환해줌.

FilterRes = FPres .* rampFilter;

% Step 3. (2) * (3) 후 inverse fourier transform.
% 이 때, (2) * (3) 해주는 방향은 detector방향.

FilterRes = ifft(FilterRes, pad_size, 2);
% Inverse Fourier Transform으로 길어진 길이는, 중앙을 기준으로 기존 sinogram과 같은 크기로 자른다.

% 자르기
FilterRes = FilterRes(:, pad_size/2+1:pad_size/2+dt);
% FilterRes = FilterRes(:, dt:2*dt+1);
% FilterRes = FilterRes(:, size(FilterRes,2)/2, size(FilterRes,2)/2 + dt);
FilterRes = arrayfun(@real, FilterRes);
subplot(1, 3, 2);
imshow(FilterRes, [0,0.005]);
title('Filtered Sinogram')

% Back Projection 구현------------------------
% FilterRes_moved = (view, dt)
% view -> (view-1) * angch
% dt -> -dt/2 ~ dt/2 -> -dt*delta/2 ~ dt*delta/2
% 함수 변환함. 이렇게 된다면 delta 때문에 보간법 필요.

px_size = 1;
FBPres = zeros(512, 512);   % 고정. 바꾸지 말 것!
[X, Y] = meshgrid(1:px_size:512, 1:px_size:512);
for k = 1:view
    t = (k-1)*angch;
    T = (X+(px_size/2)-256) * cos(t) + (Y+(px_size/2)-256 + b-SCD) * sin(t);
    kth_row = FilterRes(k,:);
    vq = interp1(1:size(kth_row, 2), kth_row, T(:)/delta + (dt-1)/2 + 1, 'linear', 0);
    %vq = interp1(1:size(kth_row, 2), kth_row, )
    vq(isnan(vq)) = 0; % NaN 값은 0.
    FBPres = FBPres + reshape(vq, size(X));
end
FBPres = (2*pi/view) * FBPres;
subplot(1, 3, 3);  
imshow(FBPres);

function [res] = FP(x, y, r, q1, q2)  %Forward Projection implementation
    C = [x;y];
    d = abs(det([q2-q1, C-q1]))/norm(q2-q1);
    if d >= r
        res = 0;
    else
        res = 2 * sqrt(r^2 - d^2);  % Find length of object at the point
    end
end

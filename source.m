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
%q1xval = [];
%q2xval= [];
%q1yval = [];
%q2yval= [];
dsxArr = [-(dt*delta/2):delta:(dt*delta/2)];  % x values of Detector & Source in Cart esian coordinate
angch = 2*pi/(view-1);
FPres = zeros(view, SDD); % Store Forward Projection Result.

cnt=0;

%Forward Projection
%Step 1. Find projection from each angle. 
for t = 0:angch:2*pi % rotate angle
    A = [cos(t) sin(t) ; -sin(t) cos(t)];
    cnt = cnt+1;
    for i = 1:dt   % #detector
        q1 = A*[dsxArr(i); SCD];
        q2 = A*[dsxArr(i); SDD-SCD];
        FPres(cnt, i) = IN*FP(a,b,r,q1,q2);
    end
end
sinogram_size = size(FPres);
subplot(1, 3, 1);   
heatmap(FPres);
title('Original Sinogram');

pad_size = 2 ^ ceil(log2(2*dt)); 
FPres = fft(FPres, pad_size, 1);

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
rampFilter = fft(rampFilter); %이 fft(X,n)을 하면 푸리에 변환 시, 함수 안에서 자체적으로 X를 n만큼 늘린 값으로 푸리에변환해줌.

FilterRes = zeros(size(FPres));
for i=1:size(FPres,2)
    for j=1:size(FPres,1)
       FilterRes(j,i) = FPres(j,i) * rampFilter(j);
    end
end

% Step 3. (2) * (3) 후 inverse fourier transform.
% 이 때, (2) * (3) 해주는 방향은 detector방향.

FilterRes = ifft(FilterRes);
% Inverse Fourier Transform으로 길어진 길이는, 중앙을 기준으로 기존 sinogram과 같은 크기로 자른다.
FilterRes = arrayfun(@real, FilterRes);
% 평행이동
FilterRes_moved = zeros(size(FilterRes));
for j=1:size(FilterRes,2)
    for i=1:size(FilterRes,1)
        FilterRes_moved(mod(i+size(FilterRes,1)/2, size(FilterRes,1)) + 1, j) = FilterRes(i, j);
    end
end
subplot(1, 3, 2);  
heatmap(FilterRes); 
title('Filtered Sinogram')

% 자르기
FilterRes_moved = FilterRes_moved(1:view, :);

subplot(1,3,3)
heatmap(FilterRes_moved);
title('Filtered Projection(Move & Cut)')



%  Back Projection 코드 ------------------------------
FBPres = zeros( ceil(dt * delta), SDD);   %가중치를 칠해줄 matrix. (x,y)이며, lower_bound 값이 본인의 좌표.)

dsx = [-(dt*delta/2):delta/5:(dt*delta/2)]; %보간법으로 만든 새 detector의 x좌표들.
dsx = dsx + (dt*delta/2);
for i = 1:view %각 view에 대하여 값 더해줄 것임.
    x = 1:dt;
    v = FilterRes(i,x);
    xq = 1:0.2:dt;
    vq = interp1(x, v, xq); %각 view에 대한 g(t, theta) interpolation
    t = (i-1)*angch;
    m = tan(t);
    A = [cos(t) sin(t); -sin(t) cos(t)];
    for j = 1:size(1:0.2:dt)  % 보간된 값에 대하여, detector가 만나는 모든 격자에 가중치를 칠해줌.
        % detector의 위치: xq
        % detector에 주어진 가중치: vq
        % 1. detector 직선함수 구하기
        % 2. 직선 함수가 지나는 index 찾아서 그 곳에 vq 더해주기
        % 개멍청한 방법이긴 한데 n^2 함수로 모든 index 돌면서 되나 안되나 찾을 수도 있음
        
        % 1. detector 직선함수 구하는 법:
        newdot = A*[dsx(j); 0];  % 회전한 detector가 지나는 하나의 점. x절편을 회전한 것이다.
        % y_dt = m * (x_range-q(1)) + q(2);   detector의 직선함수.
        cnt = 0;

        for p = 1:dt*delta
            for q = 1:SDD
                par = m * (p+0.5-newdot(1,1)) + newdot(2,1);
                if par >= q &&  par < q+1 
                    cnt = cnt+1;
                end
            end
        end

        for p = 1:dt*delta
            for q = 1:SDD
                par = m * (p+0.5-newdot(1,1)) + newdot(2,1);
                if par >= q &&  par < q+1 
                    FBPres(p,q) = FBPres(p,q) + (vq(j));
                end
            end
        end

    end
end

FBPres = (2*pi / view) * FBPres; 

function [res] = FP(x, y, r, q1, q2)  %Forward Projection implementation
    C = [x;y];
    d = abs(det([q2-q1, C-q1]))/norm(q2-q1);
    if d >= r
        res = 0;
    else
        res = 2 * sqrt(r^2 - d^2);  % Find length of object at the point
    end
end

%create the setup , click a picture and then load it 

I = imread('Indoor.jpeg');
I = rgb2gray(I);
size(I);
imshow(I);

% Select some corner points with the corner function in matlab :

%{
Minimum accepted quality of corners, specified as the comma-separated pair consisting of 
'QualityLevel' and a numeric scalar in the range (0, 1). For a quality level Q, the 
toolbox rejects candidate corners with corner metric values less than Q * max(corner metric). 
Use larger values of Q to remove erroneous corners.

Corner detection method, specified as 'Harris' for the Harris corner detector,
or 'MinimumEigenvalue' for Shi & Tomasi's minimum eigenvalue method.
%}

camCordinate = corner(I, 'MinimumEigenvalue', 100, 'QualityLevel', 0.4);
%remove some random 8 row index and keep 32 random points.

%lets display the corner points on the image using the hold on functionality 
imshow(I);
hold on
plot(camCordinate(:,1), camCordinate(:,2), 'r*');

%Lets save some coordinates for testing purpose and project them  
testcamCordinate = camCordinate([79,11,96,25,10], :);
imshow(I);
hold on
plot(testcamCordinate(:,1), testcamCordinate(:,2), 'r*');

%coresponding true world coordinates for the testing points:
trueworldCoordinate = [0*23 0*23 23*2;
    0*23 4*23 1*23;
    0*23 1*23 1*23;
    0*23 0*23 1*23;
    0*23 5*23 0*23;
    ];

% Now lets select the corners using which we will calculate the parameters
camCordinate = camCordinate([1,92,67,20,30,97,57,70,98,29,12,79,11,96,25,10,13,69,60,49,44,63,58,93,85,89,55,88,90,53,91,66], :);

%Plotting the Selected Corner Co-ordinates :
imshow(I);
hold on
plot(camCordinate(:,1), camCordinate(:,2), 'r*');

% Corresponding world cordinates of the points selected above is calculated
% manually and stored in the matrix .

worldCoordinate = [0 23 23*7;
    0 4*23 23*7;
    0 23*2 23*7;
    0 4*23 6*23;
    0 2*23 6*23;
    0 4*23 5*23;
    0 2*23 5*23;
    0 3*23 3*23;
    0 1*23 3*23;
    0 5*23 2*23;
    0 3*23 2*23;
    0 0 2*23;
    0 4*23 1*23;
    0 1*23 1*23;
    0 0 1*23;
    0 5*23 0;
    0 3*23 0;
    0 2*23 0;
    0 0 23*8;
    23 0 8*23;
    3*23 0 8*23;
    2*23 0 7*23;
    4*23 0 7*23;
    6*23 0 7*23;
    1*23 0 6*23;
    5*23 0 6*23;
    3*23 0 3*23;
    4*23 0 2*23;
    1*23 0 1*23;
    3*23 0 1*23;
    5*23 0 1*23;
    4*23 0 0;
    ];

% Now lets create the the homogeneous matrix

% create a matrix of size : (for n points) 2n*12
rows = 32; 
P(2 * rows,12) = 0;

%lets populate the values in the matrix: 
j = 1;

for i = 1:1:64
    if mod(i,2) ~= 0
        P(i,1:3) = worldCoordinate(j,1:3);
        P(i,4) = 1;
        P(i,9:12) = P(i,1:4)*-1*camCordinate(j,1);
    else
        P(i,5:7) = worldCoordinate(j,1:3);
        P(i,8) = 1;
        P(i,9:12) = P(i,5:8)*-1*camCordinate(j,2);
        j = j + 1 ;
        
    end
    
end


% singular value decomposition to find 'm' which is the singular column of
%'V' with minimum values.
%{
1. To find the solution of Ax = 0 
2. We need to find the non-trivial null vector of A.
3. ?A? can have up to 12 eigenvalues.
   a) Case 1: If rank(A) is 12, its nullity is zero. There is no non-trivial null vector of A.  
   b) Case 2: If rank(A) is 11, it will have exactly one zero eigenvalue and the corresponding 
      eigenvector will be the solution of Ax = 0 .
   c) Case 3: If rank(A) < 11, there are infinite solutions to Check if data is degenerate. 
      Recalibrate.


%}
U = transpose(P)*P;

%{
[V,D] = eig(A) returns diagonal matrix D of eigenvalues and
matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
%}
[V,D] = eig(U);

%{
Always check the smallest eigenvalue and/or the ratio
  between the largest and smallest values to estimate noise in the data.
%}

%Now lets create the M matrix using the first col of the Matrix V 
% First col is the coloumn against teh max value eighen Value 

M = zeros(3,4);
j = 1;
for i = 1:1:3
    M(i,1:4) = V(j:j+3,1);
    j = j+4;
end

A1(1:3) = M(1,1:3);
A2(1:3) = M(2,1:3);
A3(1:3) = M(3,1:3);

B = M(:,4);

% Now the 3rd row in the M matrix ix the rotation matrix R3 

norm_A3 = norm(A3);
R3 = (1/norm_A3) * A3; % rotation matrix1

% Rotation matrix R2
y_0_unscaled = dot(A2,transpose(R3));
R2_betaBySinTheta = A2 - y_0_unscaled * R3; 
norm_R2_betaBySinTheta = norm(R2_betaBySinTheta);
R2 = (1/norm_R2_betaBySinTheta) * R2_betaBySinTheta;

% temp = A2 - ( dot(A2,transpose(R3)) ) * R3;
% R2 = temp/norm(temp) ;

% Rotation Matrix R1 : 
x_0_unscaled = dot(A1, transpose(R3)); 
alpha_cotTheta_unscaled = dot(A1, transpose(R2));
alpha_R1 = A1 - (alpha_cotTheta_unscaled * R2) - (x_0_unscaled * R3);
norm_alpha_R1 = norm(alpha_R1);
R1 = (1/norm_alpha_R1) * alpha_R1; 

% Normalizing values with respect to norm_A3 since norm_A3 = 1 on matrix K.
% These are Intrinsic parameters.

alpha = norm_alpha_R1 / norm_A3;
skew = alpha_cotTheta_unscaled / norm_A3;
x_0 = x_0_unscaled / norm_A3;
beta_by_sinTheta = norm_R2_betaBySinTheta / norm_A3;
y_0 = y_0_unscaled / norm_A3;
norm_A3 = norm_A3 / norm_A3;

K =[alpha skew x_0;
    0 beta_by_sinTheta y_0;
    0 0 norm_A3;];  %Intrinsic Parameter Matrix.


theta = acos (-1*cross(A1,A3)*transpose(cross(A2,A3))/(norm(cross(A1,A3)*norm(cross(A2,A3))))); %theta comes out as 1.5 radians which is app. 86 degrees.
theta_deg = 180*7*theta/22 ;
sinTheta = sin(theta);
cosTheta = cos(theta);

%Extrinsic Parameter : Rotation Matrix 
R = [R1 R2 R3];
R = transpose(reshape(R,3,3));

t = inv(K) * B; %Translation matrix

% ------------ Calculation of reprojection Error ----------------

no_of_test_points = 5 ;
%converting the true world co-ordinate to the Homogeneous system 
newcol = ones(5,1);
final = [trueworldCoordinate newcol];
semfinal = transpose(final);

%reconstruction of the pixel values :

recreatedpixel = M * semfinal;
finalrec = transpose(recreatedpixel);
finalrec = finalrec * -1; 

for i = 1:1:no_of_test_points
    finalrec(i,1) = finalrec(i,1)/finalrec(i,3);
    finalrec(i,2) = finalrec(i,2)/finalrec(i,3);
    finalrec(i,3) = finalrec(i,3)/finalrec(i,3);
end

% finding the euclidean reprojection Error 
sum = 0 ;

for i = 1:1:no_of_test_points
    x_diff = finalrec(i,1)- testcamCordinate(i,1);
    square_x_diff = x_diff ^ 2;
    y_diff = finalrec(i,2)- testcamCordinate(i,2);
    square_y_diff = y_diff ^ 2;
    total_diff = square_x_diff + square_y_diff;
    euclid_dis = total_diff ^(0.5);
    sum = sum + euclid_dis;
end

% Average reprojection error :
average_reproj_error = sum /no_of_test_points;
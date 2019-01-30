
%create the setup , click a picture and then load it 

I = imread('Img.jpg');
I = rgb2gray(I);
size(I);
%imshow(I);

% Select some corner points with the corner function in matlab :

%{
Minimum accepted quality of corners, specified as the comma-separated pair consisting of 
'QualityLevel' and a numeric scalar in the range (0, 1). For a quality level Q, the 
toolbox rejects candidate corners with corner metric values less than Q * max(corner metric). 
Use larger values of Q to remove erroneous corners.

Corner detection method, specified as 'Harris' for the Harris corner detector,
or 'MinimumEigenvalue' for Shi & Tomasi's minimum eigenvalue method.
%}

camCordinate = corner(I, 'MinimumEigenvalue', 40, 'QualityLevel', 0.4);
%remove some random 8 row index and keep 32 random points.

camCordinate([5 25 26 27 30 34 36 38], : ) = [];
%lets display the corner points on the image using the hold on functionality 
imshow(I);
hold on
plot(camCordinate(:,1), camCordinate(:,2), 'r*');


% Corresponding world cordinates of the points selected above is calculated
% manually and stored in the matrix .

worldCoordinate = [0 84 56;
168 0 112;
28 0 252;
0 28 140;
168 0 84;
112 0 140;
0 0 28;
140 0 168;
168 0 168;
84 0 28;
112 0 168;
196 0 252;
168 0 196;
196 0 224;
112 0 112;
56 0 196;
112 0 252;
0 0 196;
28 0 168;
196 0 112;
84 0 224;
0 56 252;
140 0 112;
112 0 196;
168 0 140
140 0 84;
140 0 224; 
84 0 112;
0 112 28;
112 0 224;
140 0 252;
28 0 84];

% Noe lets create the the homogeneous matrix

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
sinTheta = sin(theta);
cosTheta = cos(theta);

%Extrinsic Parameter : Rotation Matrix 
R = [R1 R2 R3];
R = transpose(reshape(R,3,3));

t = inv(K) * B; %Translation matrix
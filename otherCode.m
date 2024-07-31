% Specify the filename
filename = 'test1.txt';

% Load the data from the file
data = load(filename);

% Extract x and y values into separate arrays
x_test1 = data(:, 1);  % Assuming the x values are in the first column
y_test1 = data(:, 2);  % Assuming the y values are in the second column

% Prompt text.
prompt = "Enter desired polynomial degree: ";

% Takes the desired degree value from user input. 
degree = input(prompt);

%set default arrays to calculate the sum 
disp("Following value are for test1: ")
x_square = []; 
x_cube = []; 
x_fourth = [];
x_fifth = []; 
x_sixth = []; 
x_y = [];
x2_y = [];
x3_y = []; 
y_dif = []; 
y_error = [];
sum_y = sum(y_test1);
y_mean = sum_y/length(y_test1);

%filling in those arrays
for i = 1: length(x_test1)
    x_square (i) = x_test1(i)*x_test1(i);
    x_cube (i) = x_test1(i)*x_test1(i)*x_test1(i);
    x_fourth (i) = x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i);
    x_fifth (i) = x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i);
    x_sixth (i) = x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i)*x_test1(i);
    x_y(i) = x_test1(i)*y_test1(i);
    x2_y(i) = x_test1(i)*x_test1(i)*y_test1(i);
    x3_y(i) = x_test1(i)*x_test1(i)*x_test1(i)*y_test1(i);
    y_dif(i) = (y_test1(i)-y_mean)^2;    
end

%calculating all the summations needed for the calculations
sum_x = sum(x_test1);
sum_x_square = sum(x_square); 
sum_x_cube = sum(x_cube); 
sum_x_fourth = sum(x_fourth); 
sum_x_fifth = sum(x_fifth);
sum_x_sixth = sum(x_sixth);
sum_x_y = sum(x_y);
sum_x2_y = sum (x2_y);
sum_x3_y = sum (x3_y); 
St = sum (y_dif); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if degree is one, fill in in its respective matrix to solve for the
%coefficients 
if degree == 1
    A = [length(x_test1), sum_x; sum_x, sum_x_square];

    B = [sum_y; sum_x_y];

%if degree is two fill in in its respective matrix to solve for the
%coefficients 
elseif degree ==2
    
    A = [length(x_test1), sum_x, sum_x_square; 
        sum_x, sum_x_square, sum_x_cube;
        sum_x_square, sum_x_cube, sum_x_fourth];

    B = [sum_y; sum_x_y; sum_x2_y]; 
 
%if degree is three fill in in its respective matrix to solve for the
%coefficients 
else 
    A = [length(x_test1), sum_x, sum_x_square, sum_x_cube; 
        sum_x, sum_x_square, sum_x_cube, sum_x_fourth;
        sum_x_square, sum_x_cube, sum_x_fourth, sum_x_fifth;
        sum_x_cube, sum_x_fourth, sum_x_fifth, sum_x_sixth];

    B = [sum_y; sum_x_y; sum_x2_y; sum_x3_y]; 
end

%***** Gaussian elimination from A1*********
% Creates augmented matrix.
aug = [A, B];

%Identifies the size of the matrix that needs to be created.
[rows, columns] = size(A);

y = 2;

% Declared shift variables.
shiftdw = 1;
shiftri = 1;
shifty = 2;
shiftI = 1;

% x counts up to column size. 
for x = 1:columns
    
    %Identifies what is the highest number in a column and swaps it with
    %the current row that we are trying to subtract. (Performs the partial pivoting)
    for i = shiftI:rows
        if  aug(i, x) == max(aug(shiftI:rows, x))
            if (aug(x, x) < max(aug(shiftI:rows, x)))
                aug([x, i], :) = aug([i, x], :);
                shiftI = shiftI + 1;
            else 
                shiftI = shiftI + 1;
                break;
            end
        end
    end

    %Performs the guassian elmination. After row has been partial pivoted
    %we subtract all rows underneath by following formula. This formula
    %applies the forward elimination. 
    for y = shifty:rows
        const = aug(y,x)/aug(shiftdw, shiftri);
        newRow = aug(y, :) - const*aug(shiftdw, :);
        aug(y, :) = newRow;        
    end

    % Shifts the position to the right and and down by one.
    if (y == rows)
       shiftdw = shiftdw + 1;
       shiftri = shiftri + 1;
       shifty = shifty + 1;
    end
    
end


% Delcares variables.
z = columns;
augSize = columns;

% Created a variable Stores final answers
coefficients = zeros(columns, 1);

%Performs back substitution and finds the x values for the matrix. It does
%this by isolating for the first x value at the last row of the matrix and
%incremently moves back each row and subs in the know values to find the
%unkown in each row.
while z ~= 0

    coefficients(z) = aug(z, augSize+1) / aug(z, z);
     
    k = z - 1;

    while k ~= 0
        aug(k, augSize+1) = aug(k, augSize+1) - aug(k, z)*coefficients(z);
        aug(k, z) = 0;
        k = k - 1;
    end
    
    z = z - 1;
end

%************************************************************************
disp("The coefficients of the polynomial are:");
disp(coefficients);

%display polynomials according to the chosen functions
disp("The function will be: ");
if degree == 1
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x']);

    %calculating e^2 as the difference between the actual and the estimate
    %at every point 
    y_error = (y_test1 - coefficients(1)*ones(length(x_test1), 1) - coefficients(2)*x_test1).^2; 

elseif degree == 2
    %calculating e^2 as the difference between the actual and the estimate
    %at every point 
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x + ' num2str(coefficients(3)) 'x^2']);
    y_error = (y_test1 - coefficients(1)*ones(length(x_test1), 1) - coefficients(2)*x_test1 - coefficients(3)*(x_test1.^2)).^2; 
else
    %calculating e^2 as the difference between the actual and the estimate
    %at every point 
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x + ' num2str(coefficients(3)) 'x^2 + ' num2str(coefficients(4)) 'x^3']);
    y_error = (y_test1 - coefficients(1)*ones(length(x_test1), 1) - coefficients(2)*x_test1 - coefficients(3)*(x_test1.^2) - coefficients(4)*(x_test1.^3)).^2; 
end

%calculating r^2
r_square = (St - sum(y_error))/St;

disp("The r^2 value is: ")
disp(r_square);

% Plot the points
figure;
scatter(x_test1, y_test1, 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
hold on;  % Keep the current plot while adding the polynomial plot

% Plot the polynomial function
polynomial_values = polyval(flip(coefficients), x_test1);
plot(x_test1, polynomial_values, 'LineWidth', 2, 'Color', 'red');

% Add labels and title
title('Regression Analysis for Test1');
xlabel('x');
ylabel('y');
legend('Raw data', 'Estimated Function');
text(x_test1(end-2),y_test1(end-4), ['R^2 = ' num2str(r_square)]);
grid on;
hold off;  % Release the current plot hold

disp("*****************************************************************************************************")
disp("Following results are for test2")

% Specify the filename
filename = 'test2.txt';

% Load the data from the file
data = load(filename);

% Extract x and y values into separate arrays
x_test2 = data(:, 1);  % Assuming the x values are in the first column
y_test2 = data(:, 2);  % Assuming the y values are in the second column

 
x_square = []; 
x_cube = []; 
x_fourth = [];
x_fifth = []; 
x_sixth = []; 
x_y = [];
x2_y = [];
x3_y = []; 
y_dif = []; 
y_error = [];
sum_y = sum(y_test2);
y_mean = sum_y/length(y_test2);


for i = 1: length(x_test2)
    x_square (i) = x_test2(i)*x_test2(i);
    x_cube (i) = x_test2(i)*x_test2(i)*x_test2(i);
    x_fourth (i) = x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i);
    x_fifth (i) = x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i);
    x_sixth (i) = x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i)*x_test2(i);
    x_y(i) = x_test2(i)*y_test2(i);
    x2_y(i) = x_test2(i)*x_test2(i)*y_test2(i);
    x3_y(i) = x_test2(i)*x_test2(i)*x_test2(i)*y_test2(i);
    y_dif(i) = (y_test2(i)-y_mean)^2;    
end

sum_x = sum(x_test2);
sum_x_square = sum(x_square); 
sum_x_cube = sum(x_cube); 
sum_x_fourth = sum(x_fourth); 
sum_x_fifth = sum(x_fifth);
sum_x_sixth = sum(x_sixth);
sum_x_y = sum(x_y);
sum_x2_y = sum (x2_y);
sum_x3_y = sum (x3_y); 
St = sum (y_dif); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if degree == 1
    A = [length(x_test2), sum_x; sum_x, sum_x_square];

    B = [sum_y; sum_x_y];
elseif degree ==2
    
    A = [length(x_test2), sum_x, sum_x_square; 
        sum_x, sum_x_square, sum_x_cube;
        sum_x_square, sum_x_cube, sum_x_fourth];

    B = [sum_y; sum_x_y; sum_x2_y]; 

else 
    A = [length(x_test2), sum_x, sum_x_square, sum_x_cube; 
        sum_x, sum_x_square, sum_x_cube, sum_x_fourth;
        sum_x_square, sum_x_cube, sum_x_fourth, sum_x_fifth;
        sum_x_cube, sum_x_fourth, sum_x_fifth, sum_x_sixth];

    B = [sum_y; sum_x_y; sum_x2_y; sum_x3_y]; 
end

% Creates augmented matrix.
aug = [A, B];

%Identifies the size of the matrix that needs to be created.
[rows, columns] = size(A);

y = 2;

% Declared shift variables.
shiftdw = 1;
shiftri = 1;
shifty = 2;
shiftI = 1;

% x counts up to column size. 
for x = 1:columns
    
    %Identifies what is the highest number in a column and swaps it with
    %the current row that we are trying to subtract. (Performs the partial pivoting)
    for i = shiftI:rows
        if  aug(i, x) == max(aug(shiftI:rows, x))
            if (aug(x, x) < max(aug(shiftI:rows, x)))
                aug([x, i], :) = aug([i, x], :);
                shiftI = shiftI + 1;
            else 
                shiftI = shiftI + 1;
                break;
            end
        end
    end

    %Performs the guassian elmination. After row has been partial pivoted
    %we subtract all rows underneath by following formula. This formula
    %applies the forward elimination. 
    for y = shifty:rows
        const = aug(y,x)/aug(shiftdw, shiftri);
        newRow = aug(y, :) - const*aug(shiftdw, :);
        aug(y, :) = newRow;        
    end

    % Shifts the position to the right and and down by one.
    if (y == rows)
       shiftdw = shiftdw + 1;
       shiftri = shiftri + 1;
       shifty = shifty + 1;
    end
    
end


% Delcares variables.
z = columns;
augSize = columns;

% Created a variable Stores final answers
coefficients = zeros(columns, 1);

%Performs back substitution and finds the x values for the matrix. It does
%this by isolating for the first x value at the last row of the matrix and
%incremently moves back each row and subs in the know values to find the
%unkown in each row.
while z ~= 0

    coefficients(z) = aug(z, augSize+1) / aug(z, z);
     
    k = z - 1;

    while k ~= 0
        aug(k, augSize+1) = aug(k, augSize+1) - aug(k, z)*coefficients(z);
        aug(k, z) = 0;
        k = k - 1;
    end
    
    z = z - 1;
end


disp("The coefficients of the polynomial are:");
disp(coefficients);

disp("The function will be: ");
if degree == 1
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x']);
    y_error = (y_test2 - coefficients(1)*ones(length(x_test2), 1) - coefficients(2)*x_test2).^2; 

elseif degree == 2
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x + ' num2str(coefficients(3)) 'x^2']);
    y_error = (y_test2 - coefficients(1)*ones(length(x_test2), 1) - coefficients(2)*x_test2 - coefficients(3)*(x_test2.^2)).^2; 
else
    disp(['y = ' num2str(coefficients(1)) ' + ' num2str(coefficients(2)) 'x + ' num2str(coefficients(3)) 'x^2 + ' num2str(coefficients(4)) 'x^3']);
    y_error = (y_test2 - coefficients(1)*ones(length(x_test2), 1) - coefficients(2)*x_test2 - coefficients(3)*(x_test2.^2) - coefficients(4)*(x_test2.^3)).^2; 
end


r_square = (St - sum(y_error))/St;

disp("The r^2 value is: ")
disp(r_square);

% Plot the points
figure;
scatter(x_test2, y_test2, 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'blue');
hold on;  % Keep the current plot while adding the polynomial plot

% Plot the polynomial function
polynomial_values = polyval(flip(coefficients), x_test2);
plot(x_test2, polynomial_values, 'LineWidth', 2, 'Color', 'red');

% Add labels and title
title('Regression Analysis for Test2');
xlabel('x');
ylabel('y');
legend('Raw data', 'Estimated Function');
text(x_test2(end-2),y_test2(end-4), ['R^2 = ' num2str(r_square)]);
grid on;
hold off;  % Release the current plot hold



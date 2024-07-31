A = readmatrix("test1.txt");

% Displays the instructions to terminal
disp("Select the function to fit your data");
disp("1. Polynomial");
disp("2. Exponential");
disp("3. Saturation");

%Asks for user input
choice = input("Enter the function number (1, 2, or 3): ");

h = 0;

% Runs code twice to display test1 and test2 data
while h < 2

    if choice == 2

    i = 1;
    x = 0;
    y = 0;
    MulThenSum = 0;
    xSqrThenSum = 0;
    a = 0;
    b = 0;
    yAvg = 0; 
    sr = 0;
    st = 0;
    
    % Creates matrix
    [rows, columns] = size(A);

    % Performs summations of X and Y
    while i < rows + 1
    
        x = x + A(i, 1);
        y = y + log(A(i, 2));

        MulThenSum = MulThenSum + log(A(i, 2))*A(i, 1);
        xSqrThenSum = xSqrThenSum + A(i, 1)*A(i, 1);
        yAvg = yAvg + A(i, 2);
        i = i + 1;
    end

    % Calculates constant a and b
    a = (rows*MulThenSum - x*y) / ((rows*xSqrThenSum) - (x*x));
    b = (y / rows) - a*(x / rows);

    B = exp(b);

    % Finds mean of y values
    yAvg = yAvg / rows;

    i = 1;

    % Calculates the coffiecent of determination Sr and St
    while i < rows + 1

        sr = sr + (A(i, 2) - B*exp(a*A(i, 1)))^2;
        st = st + (A(i, 2) - yAvg)^2;
        i = i + 1;
    end

    % Finds R^2
    rSqr = (st - sr) / st;
    
    xAxis = A(:, 1);
    yAxis = A(:, 2);

    % Caculates the exponent function for graph visual
    yFunc = B*(exp(a*xAxis));

    % Displays answers and plots data points and the exponent graph for
    % both data tables
    if (h == 1)
        
        disp("test2:");
        disp("Exponential, y = " + B + " * exp(" + a + " * x), R^2 = " + rSqr);

        figure;
        plot(xAxis, yAxis, 'ks');

        hold on

        plot(xAxis, yFunc, 'r');

        title('Regression Analysis (Exponential) for test2.txt');
        xlabel('x');
        ylabel('y');
        text(A(1, 1), A(rows, 2)*0.90, "Exponential, y = " + B + " * exp(" + a + " * x), R^2 = " + rSqr);
        legend('Raw Data', 'Estimated Function');
        hold off;  % Release the current plot hold
    else

        disp("test1:");
        disp("Exponential, y = " + B + " * exp(" + a + " * x), R^2 = " + rSqr);

        plot(xAxis, yAxis, 'ks');

        hold on

        plot(xAxis, yFunc, 'r');

        title('Regression Analysis (Exponential) for test1.txt');
        xlabel('x');
        ylabel('y');
        text(A(1, 1), A(rows, 2)*0.90, "Exponential, y = " + B + " * exp(" + a + " * x), R^2 = " + rSqr);
        legend('Raw Data', 'Estimated Function');
        hold off;  % Release the current plot hold
    end

    elseif choice == 3

    i = 1;
    x = 0;
    y = 0;
    MulThenSum = 0;
    xSqrThenSum = 0;
    a = 0;
    b = 0;
    sr = 0;
    st = 0;
    rSqr = 0;
    yAvg = 0;

    % Creates matrix
    [rows, columns] = size(A);

    % Performs summations of X and Y
    while i < rows + 1

        x = x + 1 / A(i, 1);
        y = y + 1 / A(i, 2);

        MulThenSum = MulThenSum + (1 / A(i, 2))*(1 / A(i, 1));
        xSqrThenSum = xSqrThenSum + (1 / A(i, 1))*(1 / A(i, 1));
        yAvg = yAvg + A(i, 2);
        i = i + 1;
    end

    % Calculates constant a and b
    a = (rows*MulThenSum - x*y) / ((rows*xSqrThenSum) - (x*x));
    b = (y / rows) - a*(x / rows);

    % Finds mean of y values
    yAvg = yAvg / rows;

    SatA = 1 / b;
    SatB = a*SatA;

    i = 1;

    % Calculates the coffiecent of determination Sr and St
    while i < rows + 1

        sr = sr + (A(i, 2) - ((SatA*A(i, 1)) / (SatB + A(i, 1))) )^2;
        st = st + (A(i, 2) - yAvg)^2;
        i = i + 1;
    end

    % Finds R^2
    rSqr = 1 - (sr / st);

    xAxis = A(:, 1);
    yAxis = A(:, 2);

    % Caculates the saturation function for graph visual
    satFunc = (SatA*xAxis)./(SatB + xAxis);

    % Displays answers and plots data points and the saturation graph for
    % both data tables
    if h == 1

        disp("test2:");
        disp("Saturation, y = (" + SatA + " x) / (" + SatB + " + x), R^2 = " + rSqr);

        figure;

        plot(xAxis, yAxis, 'ks');

        hold on

        plot(xAxis, satFunc, 'r');

        title('Regression Analysis (Saturation) for test2.txt');
        xlabel('x');
        ylabel('y');
        text(A(1, 1), A(rows, 2)*0.90, "Saturation, y = (" + SatA + " x) / (" + SatB + " + x), R^2 = " + rSqr);
        legend('Raw Data', 'Estimated Function');
        hold off;  % Release the current plot hold
    else
    
        disp("test1:");
        disp("Saturation, y = (" + SatA + " x) / (" + SatB + " + x), R^2 = " + rSqr);

        plot(xAxis, yAxis, 'ks');
        
        hold on

        plot(xAxis, satFunc, 'r');

        title('Regression Analysis (Saturation) for test1.txt');
        xlabel('x');
        ylabel('y');
        text(A(1, 1), A(rows, 2)*0.85, "Saturation, y = (" + SatA + " x) / (" + SatB + " + x), R^2 = " + rSqr);
        legend('Raw Data', 'Estimated Function');
        hold off;  % Release the current plot hold
    end

    else

    if h == 0

    % Specify the filename
    % filename = 'test1.txt';

    % Load the data from the file
    data = readmatrix("test1.txt");
    
    % Extract x and y values into separate arrays
    x_test1 = data(:, 1);  % Assuming the x values are in the first column
    y_test1 = data(:, 2);  % Assuming the y values are in the second column

        % Prompt text.
    prompt = "Enter desired polynomial degree: ";

    % Takes the desired degree value from user input. 
    degree = input(prompt);
    end 



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

    if h == 1
    % Add labels and title
    title('Regression Analysis for Test2');
    xlabel('x');
    ylabel('y');
    legend('Raw data', 'Estimated Function');
    text(x_test1(end-2),y_test1(end-4), ['R^2 = ' num2str(r_square)]);
    grid on;
    hold off;  % Release the current plot hold
    else
    % Add labels and title
    title('Regression Analysis for Test1');
    xlabel('x');
    ylabel('y');
    legend('Raw data', 'Estimated Function');
    text(x_test1(end-2),y_test1(end-4), ['R^2 = ' num2str(r_square)]);
    grid on;
    hold off;  % Release the current plot hold
    end
    end

h = h + 1;

if h == 1 

    % Loads the second set of data points
    A = readmatrix("test2.txt");
    data = readmatrix("test2.txt");

    % Extract x and y values into separate arrays
    x_test1 = data(:, 1);  % Assuming the x values are in the first column
    y_test1 = data(:, 2);  % Assuming the y values are in the second column
end

end



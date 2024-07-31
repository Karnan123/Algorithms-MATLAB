% Loads A & B text files.
A = load('A.txt');
B = load('B.txt');

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

    %Displays the current matrix. 
    disp(aug);
    
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

%Displays matrix.
disp(aug);

% Delcares variables.
z = columns;
augSize = columns;

% Created a variable Stores final answers
arrayX = zeros(columns, 1);

%Performs back substitution and finds the x values for the matrix. It does
%this by isolating for the first x value at the last row of the matrix and
%incremently moves back each row and subs in the know values to find the
%unkown in each row.
while z ~= 0

    arrayX(z) = aug(z, augSize+1) / aug(z, z);
     
    k = z - 1;

    while k ~= 0
        aug(k, augSize+1) = aug(k, augSize+1) - aug(k, z)*arrayX(z);
        aug(k, z) = 0;
        k = k - 1;
    end
    
    %Displays list of solutions to the x values in the matrix. 
    disp(arrayX(z))
    
    z = z - 1;
end

% Prints the final matrix.
disp(aug);

%Prints out final matrix x values. Rounds all answers to 5 digits.
arrayX = round (arrayX*1000, 5);

disp(arrayX);

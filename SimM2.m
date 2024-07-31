% Given coefficients from example 9 of "Roots of nonlinear equations".
A = 3.9083*10^-3;
B = -5.775*10^-7;
C = -4.183*10^-12;
T = 300;

% Resistance constant.
R1 = 100;

% Prompt text.
prompt = "Enter a resistance value: ";

% Intial values for iteration counter and approximate error percentage.
count = 0;
er = 100;

% Takes the resistance value from user input. 
R = input(prompt);

% Prompt text.
prompt = "Enter a relative approximate error percentage: ";

% Takes the approximate error percentage from user input.
userEr = input(prompt);

% Uses Hot RTD formula if user resistance input is greater or equal to 100
% else uses Cold RTD formula.
if R >= 100 

    % Compares the current error percentage with the desired error percentage
    % by the user.
    while er >= userEr

        % Holds the previous iteration of temperature.
        OldT = T;

        % Performs the Fixed point iteration by isolating for temperature.
        T = (R - R1 - R1*B*(T^2))/(R1*A);

        % Calculates the absolute relative approximate error percentage 
        % using the current and previous iteration. 
        er = (abs(T - OldT)/abs(T))*100;

        % Tracks how many iterations were performed.
        count = count + 1;
    end

    % Given initial temperature.
    T2 = 300;

    % Intial values for iteration counter and approximate error percentage.
    er2 = 100;
    count2 = 0;

    % Compares the current error percentage with the desired error percentage
    % by the user.
    while er2 >= userEr
        
        % Holds the previous iteration of temperature.
        OldT2 = T2;

        % Performs the Newton Rasphson method by subtracting the previous 
        % iteration of temperature from the previous iteration of the 
        % temperature function over its derivative.
        T2 = T2 - (B*(T2^2) + A*T2 + (R1 - R)/100)/(2*B*T2 + A);
    
        % Calculates the absolute relative approximate error percentage 
        % using the current and previous iteration.
        er2 = (abs(T2 - OldT2)/abs(T2))*100;

        % Tracks how many iterations were performed.
        count2 = count2 + 1;
    end

    %declare starting variables and the initial range
    Tbegin = 0;
    Tend = 850; 
    Told=0; 
    Tnew = 0; 
    biCount = 0;
    er3 = 100;

    %while loop keeps running to make sure we reach the desired approximate
    %absolute relative error for the root
    while er3 >= userEr

    %Respective value of functions at the end of the intervals
    Fbegin =  R - R1*(1 + A*Tbegin + B*Tbegin*Tbegin);
    Fend = R - R1*(1 + A*Tend + B*Tend*Tend);

    %check if the F(tnew) * F(end) is negative to make sure the root is
    %between the interval
    if Fbegin*Fend < 0 

        %set new t which is the average of the interval
        Tnew = (Tend+Tbegin)/2;
        
        %value of the function at Tnew
        Fnew = R - R1*(1 + A*Tnew + B*Tnew*Tnew);

        %calculating absolute respective error
        er3 = 100*abs((Tnew - Told)/Tnew);

        %set new interval by checking if F(Tnew)* [F(Tbegin) OR F(Tend)] is
        %negative
        if Fbegin*Fnew < 0 
            Tend = Tnew; 
        else  
            Tbegin = Tnew;
        end
    end 

    %update counter
    biCount = biCount+1; 

    %set previous value to current value for next iteration
    Told=Tnew;

    end

else

    while er >= userEr

        % Holds the previous iteration of temperature.
        OldT = T;

        % Performs the Fixed point iteration by isolating for temperature.
        T = (R - R1 - R1*B*(T^2) - R1*C*(T-100)*(T^3))/(R1*A);

        % Calculates the absolute relative approximate error percentage 
        % using the current and previous iteration. 
        er = (abs(T - OldT)/abs(T))*100;

        % Tracks how many iterations were performed.
        count = count + 1;
    end

    % Given initial temperature.
    T2 = 300;

    % Intial values for iteration counter and approximate error percentage.
    er2 = 100;
    count2 = 0;

    % Compares the current error percentage with the desired error percentage
    % by the user.
    while er2 >= userEr
        
        % Holds the previous iteration of temperature.
        OldT2 = T2;

        % Performs the Newton Rasphson method by subtracting the previous 
        % iteration of temperature from the previous iteration of the 
        % temperature function over its derivative.
        T2 = T2 - (B*(T2^2) + A*T2 + C*(T2-100)*(T2^3) + (R1 - R)/100)/(A + 2*B*T2 + 4*C*(T2^3) - 300*C*(T2^2));
    
        % Calculates the absolute relative approximate error percentage 
        % using the current and previous iteration.
        er2 = (abs(T2 - OldT2)/abs(T2))*100;

        % Tracks how many iterations were performed.
        count2 = count2 + 1;
    end

        %declare starting variables and the initial range
    Tbegin = -200;
    Tend = -1; 
    Told=0; 
    Tnew = 0; 
    biCount = 0;
    er3 = 100;

    %while loop keeps running to make sure we reach the desired approximate
    %absolute relative error for the root
    while er3 >= userEr

        %Respective value of functions at the end of the intervals
        Fbegin =  R - R1*(1 + A*Tbegin + B*Tbegin*Tbegin + C*(Tbegin-100)*Tbegin*Tbegin*Tbegin);
        Fend = R - R1*(1 + A*Tend + B*Tend*Tend + C*(Tend-100)*Tend*Tend*Tend);

        %check if the F(tnew) * F(end) is negative to make sure the root is
        %between the interval
        if Fbegin*Fend < 0 
            %set new t which is the average of the interval
            Tnew = (Tend+Tbegin)/2;

            %value of the function at Tnew
            Fnew = R - R1*(1 + A*Tnew + B*Tnew*Tnew +C*(Tnew-100)*Tnew*Tnew*Tnew);

            %calculating absolute respective error
            er3 = 100*abs((Tnew - Told)/Tnew);
            %set new interval by checking if F(Tnew)* [F(Tbegin) OR F(Tend)] is
            %negative
            if Fbegin*Fnew < 0 
                Tend = Tnew; 
            else  
                Tbegin = Tnew;
            end
        end 
        %update counter
        biCount = biCount+1; 

        %set previous value to current value for next iteration
        Told=Tnew;
    end
end

% Prints out temperature to terminal.
disp("The temperature obtained by bisection is " + Tnew + "C");
disp("The temperature obtained by fixed point is " + T + "C");
disp("The temperature obtained by NR is " + T2 + "C");

% Prints out number of iterations to terminal.
fprintf('\n');
disp("The number of required iterations for bisection is " + biCount);
disp("The number of required iterations for fixed point is " + count);
disp("The number of required iterations for Newton-Raphson is " + count2);

% Prints out absolute relative approximate error to terminal.
fprintf('\n');
disp("The absolute relative approximate error % for bisection is " + er3);
disp("The absolute relative approximate error % for fixed point is " + er);
disp("The absolute relative approximate error % for NR is " + er2);

%-- .- -.. . / -... -.-- / -.-- --- --. . ... ....
%Date:10_11_2023
clc;clear all;close all;
N=300;    %number of time steps;
nu=2; ny=2; ne=0;   %As per the problem being solved
deg=3;              %Nonlinear degree 3
U=2*rand(N, 1)-1;    %acknowledged that U mean is not zero out of this way
noise=0.04*randn(N, 1);
Y=[0+noise(1); 0+noise(2)];    %Y vector initialisez with the value of y(1) and y(0)
%Data generation
for i=3:N
    Y=[Y; -0.6*Y(i-1)-0.18*Y(i-1)*Y(i-2)+0.58*U(i-1)-0.24*U(i-2)+noise(i)];
end
t=1:N;
%plot(t, Y);


%Dictionary generation
[DICT, DICT_HEADER]=dict_element({}, Y, U, noise, 1, deg, nu, ny, ne);
fprintf("length of dict=%d",length(DICT(1,:)));
fprintf("length of Y=%d",length(Y(3:N)));
%DICT=DICT';
[Y(1:N-2) U(1:N-2) noise(1:N-2) DICT];


% Initialization
max_terms = size(DICT, 2);  % Maximum number of terms in the dictionary
FROLS_tol = 1e-6;  % Tolerance for the FROLS algorithm
max_iter = 100;    % Maximum iterations for FROLS

% Initialize variables for FROLS algorithm
selected_terms = [];  % List to store the selected regressor indices
B = zeros(max_terms, max_terms);  % Matrix to store coefficients
residuals = Y;  % Initialize residuals with output

for iter = 1:max_iter
    min_err = inf;  % Initialize minimum error
    
    % Loop through each term in the dictionary to find the best regressor
    for i = 1:max_terms
        if ~ismember(i, selected_terms)
            % Include the term and solve for coefficients
            curr_terms = [selected_terms i];
            B(curr_terms, curr_terms) = pinv(DICT(:, curr_terms)' * DICT(:, curr_terms)) * DICT(:, curr_terms)' * Y;
            
            % Calculate residuals
            y_pred = DICT(:, curr_terms) * B(curr_terms, curr_terms);
            residuals_new = Y - y_pred;
            
            % Calculate error
            err = norm(residuals_new);
            
            % Check if error is reduced
            if err < min_err
                min_err = err;
                best_term = i;
            end
        end
    end
    
    % Update selected terms and residuals
    selected_terms = [selected_terms best_term];
    residuals = residuals_new;
    
    % Check convergence
    if min_err < FROLS_tol
        break;
    end
end

% The selected_terms now contain the indices of the significant regressors
% The B matrix holds the estimated coefficients



function [dict, dict_header]=dict_element(vec, Y, U, noise, ip, deg, nu, ny, ne)
    clc;
    N=length(Y);
    max_delay=max([nu, ny, ne]);
    max_N=N-max_delay;
    fprintf("Max delay found= %d Max_N= %d\n", max_delay, max_N);
    dict_header=[];
    index=0;
    for t=1:deg
        l=length(vec);
        size_vec=size(vec)
        if l==0
            for j=ip:(ny+1+nu+1+ne+1)
                dict_header=[dict_header j]
                if j<=ny+1;
                    vec=[vec; Y(max_delay+2-j:max_N+(max_delay+1-j))'];
                elseif and(j>ny+1, j<=ny+1+nu+1);
                    vec=[vec; U(max_delay+2-(j-ny-1):max_N+(max_delay+1-(j-ny-1)))'];
                    %fprintf("j=%d\n", j);
                else and(j>ny+1+nu+1, j<=ny+1+nu+1+ne+1);
                    %fprintf("Noise j=%d\n", j);
                    vec=[vec; noise(max_delay+2-(j-ny-1-nu-1):max_N+(max_delay+1-(j-ny-1-nu-1)))'];
                end
            end
        else
            new_vec={};
            fprintf("Deg = %d Length of vec=%d and size=", deg, length(vec));
            index1=0;
            new_dict_header=[];
            for i=1:(ny+1+nu+1+ne+1)        %i corresponds to top layer
                new_subvec=[];                
                sz=size(vec{i});
                for j=i:(ny+1+nu+1+ne+1)    %j corresponds to sub layer
                    fprintf("j= %d, length(vec(j))=%d, size=\n",j, length(vec{j}));
                    sz=size(vec{j});
                    index2=0;
                    for t1=1:j-1 %Just for header purpose
                        sztemp=size(vec{t1});
                        index2=index2+sztemp(1);
                    end
                    for k=1:sz(1)       %k corresponds to element in sublayer
                        index1=index1+1;
                        fprintf("index1=%d, index2=%d, k=%d", index1, index2, k);
                        new_dict_header(index1)=(10^(t-1))*i+dict_header(index2+k)
                        if t==2
                            temp=vec{j};
                        else
                            temp=vec{j}(k,:);
                            %pause(5)
                        end
                        %fprintf("temp=");
                        %temp
                        if i<=ny+1;                        
                            new_subvec=[new_subvec; temp.*Y(max_delay+2-i:max_N+(max_delay+1-i))'];
                        elseif and(i>ny+1, i<=ny+1+nu+1);
                            new_subvec=[new_subvec; temp.*U(max_delay+2-(i-ny-1):max_N+(max_delay+1-(i-ny-1)))'];
                            %fprintf("Else ip=%d j=%d\n", ip, j);
                        else and(i>ny+1+nu+1, i<=ny+1+nu+1+ne+1);
                            %fprintf("Noise j=%d\n", j)
                            new_subvec=[new_subvec; temp.*noise(max_delay+2-(i-ny-1-nu-1):max_N+(max_delay+1-(i-ny-1-nu-1)))'];                            
                        end
                        %pause(5)
                    end           
                end
                new_vec=[new_vec [new_subvec]]
            end
            vec=new_vec;
            fprintf("Abhi dekhenge")
            dict_header=new_dict_header;
        end
    deg=deg-1;
    fprintf("We got a round here with new deg=%d and results:\n\n\n\n\n\n", deg);
    dict=cellarraytomat(vec);
    end
end

function vec=cellarraytomat(vecs)
vec=[];
for i=1:length(vecs)
    vec=[vec; vecs{i}];
end
vec=vec';
end


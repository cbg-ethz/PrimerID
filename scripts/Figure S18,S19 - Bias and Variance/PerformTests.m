warning off

x0 = [...
    0.096294; 0.145240; 0.375070; 0.231780; 0.151610;
    0.096294; 0.145240; 0.375070; 0.231780; 0.151610;
    0.096294; 0.145240; 0.375070; 0.231780; 0.151610;
    1000; 1000; 1000];

options.TolFun = 1E-15;
options.TolCon = 1E-15;
options.TolX   = 1E-15;
options.MaxFunEvals = 1E6;

objFunc = @(x) (-1)*PDF_DATA(Data_raw, Data_pID, Data_HC, x);
LB = zeros(18, 1);
UB = 1E6*ones(18, 1);


%% 1.) TEST BIASES
% fully constrained, 7 d.f.
Full_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
];
Full_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Full_x, Full_fval] = fmincon(objFunc, x0, [], [], Full_Aeq, Full_beq, LB, UB, [], options);


% partially constrained, 11 d.f.
Partial_raw_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
];
Partial_raw_beq = [1; 1; 1; 0; 0; 0; 0; 0;];
[Partial_raw_x, Partial_raw_fval] = fmincon(objFunc, x0, [], [], Partial_raw_Aeq, Partial_raw_beq, LB, UB, [], options);

Partial_pID_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
];
Partial_pID_beq = [1; 1; 1; 0; 0; 0; 0; 0;];
[Partial_pID_x, Partial_pID_fval] = fmincon(objFunc, x0, [], [], Partial_pID_Aeq, Partial_pID_beq, LB, UB, [], options);

Partial_HC_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
];
Partial_HC_beq = [1; 1; 1; 0; 0; 0; 0; 0;];
[Partial_HC_x, Partial_HC_fval] = fmincon(objFunc, x0, [], [], Partial_HC_Aeq, Partial_HC_beq, LB, UB, [], options);


% unconstrained, 15 d.f.
Uncon_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
];
Uncon_beq = [1; 1; 1;];
[Uncon_x, Uncon_fval] = fmincon(objFunc, x0, [], [], Uncon_Aeq, Uncon_beq, LB, UB, [], options);


% Inference
format longG
fprintf('\n\n');
fprintf('Full constrained: %.8f\ns: ', -Full_fval);
fprintf('%.1f  ', Full_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Full_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Full_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Full_x(11:15)');
fprintf('\n\n');



Delta_df = rank(Full_Aeq)-rank(Partial_raw_Aeq);
fprintf('\nDelta d.f.: %d\n', Delta_df);

fprintf('raw unconstrained: %.8f\n', -Partial_raw_fval);
G = -2*(-Full_fval + Partial_raw_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.1f  ', Partial_raw_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_raw_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_raw_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_raw_x(11:15)');
fprintf('\n\n');

fprintf('pID unconstrained: %.8f\n', -Partial_pID_fval);
G = -2*(-Full_fval + Partial_pID_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.1f  ', Partial_pID_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_pID_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_pID_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_pID_x(11:15)');
fprintf('\n\n');

fprintf('HC unconstrained: %.8f\n', -Partial_HC_fval);
G = -2*(-Full_fval + Partial_HC_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.1f  ', Partial_HC_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_HC_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_HC_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_HC_x(11:15)');
fprintf('\n\n');



Delta_df = rank(Full_Aeq)-rank(Uncon_Aeq);
fprintf('\nDelta d.f.: %d\n', Delta_df);
fprintf('Full unconstrained: %.8f\n', -Uncon_fval);
G = -2*(-Full_fval + Uncon_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.1f  ', Uncon_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Uncon_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Uncon_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Uncon_x(11:15)');
fprintf('\n');





%% 2.) TEST VARIANCES
% fully constrained, 5 d.f.
Full_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % s_raw == s_pID
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0;
    % s_raw == s_HC
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1;
    % s_pID == s_HC
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1;
];
Full_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Full_x, Full_fval] = fmincon(objFunc, x0, [], [], Full_Aeq, Full_beq, LB, UB, [], options);


% partially constrained, 6 d.f.
Partial_raw_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % s_pID == s_HC
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1;
];
Partial_raw_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Partial_raw_x, Partial_raw_fval] = fmincon(objFunc, x0, [], [], Partial_raw_Aeq, Partial_raw_beq, LB, UB, [], options);

Partial_pID_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % s_raw == s_HC
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1;
];
Partial_pID_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Partial_pID_x, Partial_pID_fval] = fmincon(objFunc, x0, [], [], Partial_pID_Aeq, Partial_pID_beq, LB, UB, [], options);

Partial_HC_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % s_raw == s_pID
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0;
];
Partial_HC_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Partial_HC_x, Partial_HC_fval] = fmincon(objFunc, x0, [], [], Partial_HC_Aeq, Partial_HC_beq, LB, UB, [], options);


% unconstrained, 7 d.f.
Uncon_Aeq = [...
    1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; % mu_raw sum to 1
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0; % mu_pID sum to 1
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0; % mu_HC  sum to 1
    
    % mu_raw == mu_pID
    1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0;
    
    % mu_pID == mu_HC
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0;
    
    % mu_raw == mu_HC
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0;
];
Uncon_beq = [1; 1; 1; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;];
[Uncon_x, Uncon_fval] = fmincon(objFunc, x0, [], [], Uncon_Aeq, Uncon_beq, LB, UB, [], options);


% Inference
format longG
fprintf('\n\n');
fprintf('Full constrained: %.8f\ns: ', -Full_fval);
fprintf('%.4f  ', Full_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Full_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Full_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Full_x(11:15)');
fprintf('\n\n');



Delta_df = rank(Full_Aeq)-rank(Partial_raw_Aeq);
fprintf('\nDelta d.f.: %d\n', Delta_df);

fprintf('raw unconstrained: %.8f\n', -Partial_raw_fval);
G = -2*(-Full_fval + Partial_raw_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.4f  ', Partial_raw_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_raw_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_raw_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_raw_x(11:15)');
fprintf('\n\n');

fprintf('pID unconstrained: %.8f\n', -Partial_pID_fval);
G = -2*(-Full_fval + Partial_pID_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.4f  ', Partial_pID_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_pID_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_pID_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_pID_x(11:15)');
fprintf('\n\n');

fprintf('HC unconstrained: %.8f\n', -Partial_HC_fval);
G = -2*(-Full_fval + Partial_HC_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.4f  ', Partial_HC_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Partial_HC_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Partial_HC_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Partial_HC_x(11:15)');
fprintf('\n\n');



Delta_df = rank(Full_Aeq)-rank(Uncon_Aeq);
fprintf('\nDelta d.f.: %d\n', Delta_df);
fprintf('Full unconstrained: %.8f\n', -Uncon_fval);
G = -2*(-Full_fval + Uncon_fval);
p = chi2cdf(G, Delta_df, 'upper');
fprintf('G: %.8f, p: %.8f\ns: ', G, p);
fprintf('%.4f  ', Uncon_x(16:18)');
fprintf('\nmu_raw: ');
fprintf('%.4f  ', Uncon_x(1:5)');
fprintf('\nmu_pID: ');
fprintf('%.4f  ', Uncon_x(6:10)');
fprintf('\nmu_HC : ');
fprintf('%.4f  ', Uncon_x(11:15)');
fprintf('\n');

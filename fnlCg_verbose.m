function [x,k,total_elapsed_time] = fnlCg_verbose(x0,params,options)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x)
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
%
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------

% Verbosity (written to stdout or a log file) added by BJ Anderson, CIVM,
% 27 September 2017

if exist('options','var')
    if ~isfield(options,'verbosity')
        verbosity = 0;
    else
        verbosity = options.verbosity;
        
        if ~isfield(options,'log_mode')
            log_mode = 1;
        else
            log_mode = options.log_mode;
        end
        
        if ~isfield(options,'log_file')
            log_file = '';
            log_mode = 3;
        else
            log_file = options.log_file;
        end
    end
    
    if ~isfield(options,'variable_iterations')
        variable_iterations = 0;
    else
        variable_iterations =  options.variable_iterations;
    end
    
    if (variable_iterations)
        if ~isfield(options,'convergence_limit')
            convergence_limit=0.0001;
        else
            convergence_limit=options.convergence_limit;
        end
        
        
    end
    
    if ~isfield(options,'convergence_window')
        convergence_window=10;
    else
        convergence_window=options.convergence_window;
    end
else
    verbosity=1;
    log_file='';
    log_mode=3;
    variable_iterations = 0;
    convergence_window = 10;
end

x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))

g0 = wGradient(x,params);

dx = -g0;


% BJA-2017 code
f1_vector=zeros([1 params.Itnlim]);

mean_window=zeros([1 params.Itnlim]);
convergence_metric=zeros([1 params.Itnlim]);
delta_obj = zeros([1 params.Itnlim]);
iteration_time =  zeros([1 params.Itnlim]);
max_iterations = params.Itnlim;

% iterations
if verbosity
    if (variable_iterations)
        log_msg = sprintf('Performing compressed sensing 2D slice reconstruction in VARIABLE ITERATION mode.\n\tmax iterations: %i\n\tconvergence window: %i iterations\n\tconvergence limit: %.4e\n\tconvergence metric: %s\n',max_iterations,convergence_window,convergence_limit,'Second derivative of normalized objective');
        yet_another_logger(log_msg,log_mode,log_file);
    else
        log_msg = sprintf('Performing compressed sensing 2D slice reconstruction in FIXED ITERATION mode.\n\tnumber of iterations: %i.\n',max_iterations);
        yet_another_logger(log_msg,log_mode,log_file);
    end
    
    
    log_msg = sprintf('Iteration,\t     obj,\tdelta obj,\tconv metric,\titeration time,\telapsed time.\n');
    yet_another_logger(log_msg,log_mode,log_file);
else
    time_zero = tic;
end

while(1)
    
    
    % backtracking line-search
    if verbosity
        iteration_timer = tic;
    end
    % pre-calculate values, such that it would be cheap to compute the objective
    % many times for efficient line-search
    [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
    f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, 0, params);
    t = t0;
    [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
    
    lsiter = 0;
    
    while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
        lsiter = lsiter + 1;
        t = t * beta;
        [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
    end
    
    if lsiter == maxlsiter
        disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
        return;
    end
    
    % control the number of line searches by adapting the initial step search
    if lsiter > 2
        t0 = t0 * beta;
    end
    
    if lsiter<1
        t0 = t0 / beta;
    end
    
    x = (x + t*dx);
    
    %conjugate gradient calculation
    
    g1 = wGradient(x,params);
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    g0 = g1;
    dx =  - g1 + bk* dx;
    k = k + 1;
    
    % BJA-2017 Code
    
    %RMSerr_vector(k)=RMSerr;
    
    if k==1
        f1_prime = f1;
    end
    f1_vector(k)=f1/f1_prime;
    
    if ((k >= convergence_window) && (k > 2));
        %c_window_data = RMSerr_vector((k-avgRMS_window + 1):k);
        c_window_data = f1_vector((k-convergence_window + 1):k);
        mean_window(k) = mean(c_window_data);
        %std_window(k) = std(c_window_data);
        delta_obj(k) =mean_window(k)-mean_window(k-1);
        %residual(k) = delta_mean(k)-delta_mean(k-1);%std_window(k) - abs(delta_mean_RMSerr(k));
        convergence_metric(k) = mean_window(k)-2*mean_window(k-1)+mean_window(k-2);
    end
    
    
    if (verbosity)
        iteration_time(k) = toc(iteration_timer);
        if k == 1
            elapsed_time(k) = iteration_time(k);
        else
            elapsed_time(k) = elapsed_time(k-1) + iteration_time(k);
        end

        %--------- uncomment for debug purposes ------------------------
        %log_msg = sprintf('%d   obj: %f, RMS: %f, L-S: %d\n', k,f1,RMSerr,lsiter);
        if ((k > convergence_window) && (k > 2))
            log_msg = sprintf(['%0'  num2str(numel(num2str(max_iterations)))  'i,\t%03.06f,\t%+.4e,\t%+.4e,\t%0.04f s,\t%0.04f s.\n'], k,f1,delta_obj(k), convergence_metric(k),iteration_time(k),elapsed_time(k));
        else
            log_msg = sprintf(['%0'  num2str(numel(num2str(max_iterations)))  'i,\t%03.06f,\t       %+s,\t      %+s,\t%0.04f s,\t%0.04f s.\n'], k,f1,'N/A','N/A',iteration_time(k),elapsed_time(k));
        end
        
        yet_another_logger(log_msg,log_mode,log_file);
        %---------------------------------------------------------------
    end
    
    
    
    %if (k > params.Itnlim) | (norm(dx(:)) < gradToll)
    break_flag = 0;
    if (variable_iterations)
        if ((k > convergence_window) && (k > 2) && (convergence_limit >= abs(convergence_metric(k))))
            break_flag = 1;
        end
    end
    
    if ((k >= params.Itnlim) || (break_flag))
        
        if verbosity
            total_elapsed_time = elapsed_time(k);
            if (variable_iterations)
                log_msg = sprintf('\nReconstruction complete after %i iterations; total elapsed time: %0.04f s.\nAbsolute value of the convergence metric %.4e was less than the convergence limit %.4e.\n\n',k,elapsed_time(k),convergence_metric(k),convergence_limit);
            else
                log_msg = sprintf('\nReconstruction complete after %i iterations; total elapsed time: %0.04f s.\n\n',k,elapsed_time(k));
            end
            yet_another_logger(log_msg,log_mode,log_file);
            break;
        else
            total_elapsed_time = toc(time_zero);
            break;
        end
    end
    
end


return;


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

FTXFMtx = params.FT*(params.XFM'*x);
FTXFMtdx = params.FT*(params.XFM'*dx);

if params.TVWeight
    DXFMtx = params.TV*(params.XFM'*x);
    DXFMtdx = params.TV*(params.XFM'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end





function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params);
%calculated the objective function

p = params.pNorm;

obj = FTXFMtx + t*FTXFMtdx - params.data;
obj = obj(:)'*obj(:);

if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    TV = 0;
end

if params.xfmWeight
    w = x(:) + t*dx(:);
    XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end



TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,params);
if params.xfmWeight
    gradXFM = gXFM(x,params);
end
if params.TVWeight
    gradTV = gTV(x,params);
end


grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);



function gradObj = gOBJ(x,params);
% computes the gradient of the data consistency

gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));

gradObj = 2*gradObj ;

function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);


function grad = gTV(x,params)
% compute gradient of TV operator

p = params.pNorm;

Dx = params.TV*(params.XFM'*x);

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.XFM*(params.TV'*G);







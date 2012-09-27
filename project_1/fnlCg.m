function x = fnlCg(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given 2D Fourier space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F*x - y||^2 + lambda1*TV(x)  + lambda2*f(x)...etc 
%

%-------------------------------------------------------------------------
x = x0;


% line search parameters - Dont touch..leave alone
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha; ,    beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))
g0 = wGradient(x,params);

dx = -g0;

% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	f0 = objective(x,dx, 0, params);
	t = t0;

        [f1, ERRobj, RMSerr]  =  objective(x,dx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(x,dx, t, params);
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

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d', k,f1,RMSerr,lsiter));
	%---------------------------------------------------------------
	
    %conjugate gradient calculation- Dont touch
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
    end
    
end


return;

function [res, obj, RMS] = objective(x,dx,t);
%calculated the objective function

p = params.pNorm;

%define residual

%res = obj + (TV); This was my residual for the example in class

function grad = wGradient(x,params)

%Define this function

%grad = (gradObj + params.TVWeight.*gradTV); This was my definition

function gradObj = gOBJ(x,params);
% computes the gradient of the data consistency

%Define gradient for obj function in residual
%I had defined gradient for data constraint

function grad = gTV(x,params)
% compute gradient of TV operator

%I had defined gradient for TV function

%YOU MAY WANT TO ADD MORE FUNCTIONS AND GRADIENTS






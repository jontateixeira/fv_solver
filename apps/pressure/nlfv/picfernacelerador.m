% Para uma boa convergencia tem mexer em "atol" e "rtol" a eficiencia do acelerador muito depende eles
% eu nos meus testes utilizei 1.e^-6, eu acho que você tera que fazer algum estudo parametrico para achar um 
% um adequado.

function [ x ] = picfernacelerador(x,parameter,w,s,nflagface,nflag,gamma,...
                                   weightDMP,auxface,wells,mobility,formethod)

mMax = min(10,size(x,1)); itmax = 100; DG = [];
atol = 1.e-6; rtol = 1.e-6; droptol = 1.e10; 
beta = 1; accstart = 0; res_hist = []; mAA = 0;

for iter = 0:itmax
    
    %% Interpolação das pressões na arestas (faces)
    [xinterp_new] = pressureinterp(x,nflagface,nflag,w,s,parameter,weightDMP,mobility);
                                      
    %% Calculo da matriz global
    if strcmp(formethod,'NLFVPP')==1 
        [ M_new, I_new ] = assemblematrix_NLFV(xinterp_new,parameter,wells,mobility);
    elseif strcmp(formethod,'NLFVDMP')==1 
        [ M_new, I_new ] = assemblematrixDMPSY(x,xinterp_new,gamma,parameter,weightDMP,auxface,wells,mobility);
    end
    
    gval = M_new\I_new;
    fval = gval - x;
    res_norm = norm(fval);
    res_hist = [res_hist;[iter,res_norm]];
    % Set the residual tolerance on the initial iteration.
    if iter == 0,
        tol = max(atol,rtol*res_norm);
    end
    % Test for stopping.
    if res_norm <= tol,
        break;
    end
    if mMax == 0 || iter < accstart,
        % Without acceleration, update x <- g(x) to obtain the next
        % approximate solution.
        x = gval;
    else
        % Apply Anderson acceleration.
        % Update the df vector and the DG array.
        if iter > accstart,
            df = fval-f_old;
            if mAA < mMax,
                DG = [DG gval-g_old];
            else
                DG = [DG(:,2:mAA) gval-g_old];
            end
            mAA = mAA + 1;
        end
        f_old = fval;
        g_old = gval;
        if mAA == 0
            % If mAA == 0, update x <- g(x) to obtain the next approximate solution.
            x = gval;
        else
            % If mAA > 0, solve the least-squares problem and update the
            % solution.
            if mAA == 1
                % If mAA == 1, form the initial QR decomposition.
                R(1,1) = norm(df);
                Q = R(1,1)\df;
            else
                % If mAA > 1, update the QR decomposition.
                if mAA > mMax
                    % If the column dimension of Q is mMax, delete the first column and
                    % update the decomposition.
                    [Q,R] = qrdelete(Q,R,1);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described below.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA-1); R = R(1:mAA-1,:);
                    end
                    % Explanation: If Q is not square, then qrdelete(Q,R,1) reduces the
                    % column dimension of Q by 1 and the column and row
                    % dimensions of R by 1. But if Q *is* square, then the
                    % column dimension of Q is not reduced and only the column
                    % dimension of R is reduced by one. This is to allow for
                    % MATLABâ€™s default "thick" QR decomposition, which always
                    % produces a square Q.
                end
                % Now update the QR decomposition to incorporate the new
                % column.
                for j = 1:mAA - 1
                    R(j,mAA) = Q(:,j)'*df;
                    df = df - R(j,mAA)*Q(:,j);
                end
                R(mAA,mAA) = norm(df);
                Q = [Q,R(mAA,mAA)\df];
            end
            if droptol > 0
                % Drop residuals to improve conditioning if necessary.
                condDF = cond(R);
                while condDF > droptol && mAA > 1
                    [Q,R] = qrdelete(Q,R,1);
                    DG = DG(:,2:mAA);
                    mAA = mAA - 1;
                    % The following treats the qrdelete quirk described above.
                    if size(R,1) ~= size(R,2),
                        Q = Q(:,1:mAA); R = R(1:mAA,:);
                    end
                    condDF = cond(R);
                end
            end
            % Solve the least-squares problem.
            gamma = R\(Q'*fval);
            % Update the approximate solution.
            x = gval - DG*gamma;
            % Apply damping if beta is a function handle or if beta > 0
            % (and beta ~= 1).
            if isa(beta,'function_handle'),
                x = x - (1-beta(iter))*(fval - Q*R*gamma);
            else
                if beta > 0 && beta ~= 1,
                    x = x - (1-beta)*(fval - Q*R*gamma);
                end
            end
        end
    end

end
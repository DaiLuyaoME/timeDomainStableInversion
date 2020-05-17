function [output] = stable_inversion(FF,u,L)
    FF = minreal(FF);
    if isstable(FF) == 1
        y = lsim(FF,u(1:L));
    else
        [b,a] = tfdata(FF,'v');
        
        [FA,FB,FC,FD] = tf2ss(b,a);
        [T0,lanta0] = eig(FA);
        % 纠正特征值排序
        d0 = diag(lanta0);
        d1 = abs(d0);
        [lanta,L_index] = sort(d1,'descend');
        T = T0(:,L_index);
        lanta = d0(L_index);
        
        cut = 1;
        for k = 1:length(lanta)
            if abs(lanta(k)) > 1
                cut = cut + 1;
            else
                break;
            end
        end
        d = diag(lanta);
        [TR,lantaR] = cdf2rdf(T,d); %复数特征值变为实数
        A = lantaR;
        B = pinv(TR)*FB;
        C = FC*TR;
        D = FD;
%         num = 0;
%         ET = abs(eig(TR));
%         for ii = 1:length(ET)            
%             if d1(ii) < 1e-50
%                 num = num + 1;
%             end
%             Deg = length(d1);% - num;
%             [U,S,V] = svd(TR);
%             invTR = V(:,1:Deg)*inv(S(1:Deg,1:Deg))*(U(:,1:Deg))';
%             B = invTR*FB;
%         end
        
        Aun = A(1:cut-1,1:cut-1);Bun = B(1:cut-1,:);Cun = C(:,1:cut-1);Dun = 0;
        Ast = A(cut:length(d),cut:length(d));Bst = B(cut:length(d),:);Cst = C(:,cut:length(d));Dst = D;
        N = length(u);
        xst = zeros(length(Ast),N);
        for k = 2:N
            xst(:,k) = Ast*xst(:,k-1) + Bst*u(k-1);
        end
        
        xun = zeros(length(Aun),N);
        u_L = flipud(u);
        for k = 2:N
            xun(:,k) = inv(Aun)*(xun(:,k-1)-Bun*u_L(k));
        end
        xun = fliplr(xun);
        
        y = zeros(N,1);
        for k = 1:N
            y(k) = Cst*xst(:,k) + Cun*xun(:,k) + D*u(k);
        end
    end
    output = y(1:L); %将末端发散部分截掉
end
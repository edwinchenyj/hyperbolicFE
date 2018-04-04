clear all
close all

load('eigen_decomp_bar_fix_multi_res')

lowD = sort(diag(lowD));

[highD, high_permutation_indices] = sort(diag(highD));

% lowdiag = diag(lowdiag);
% plot((lowdiag(1:10)));
plot((highD));
% title('high resolution eigenvalues of M^{-1}K')
% ylabel('eigenvalues');
% xlabel('sorted indices');

low_eig = (lowD);
high_res_low_eig = (highD);

polyfit_modes = 7;
plot_modes = 162;
low_eig(1:polyfit_modes)
high_res_low_eig(1:polyfit_modes)
% poly = polyfit([0;low_eig(1:polyfit_modes)],[0;high_res_low_eig(1:polyfit_modes)],polyfit_modes);

%%

% Minv_K_new = polyvalm(poly,low_res.M\low_res.K);

% [pV,pD] = eig(Minv_K_new);
% plowD = sort(diag(pD));
% plowD(1:polyfit_modes)



low_Minv_K = low_res.M\low_res.K;

lowdiag_unsort = lowV \ low_Minv_K * lowV;
[~, permutation_indices] = sort(diag(lowdiag_unsort));
lowV = lowV(:,permutation_indices);
lowdiag = lowV \ low_Minv_K * lowV;


lowdiag = diag(lowdiag);
plot((lowdiag(1:plot_modes)));
% plot((lowdiag));
title(['first ' num2str(plot_modes) ' low resolution eigenvalues of M^{-1}K'])
ylabel('eigenvalues');
xlabel('sorted indices');
%% Lagrange polynomial
% L = cell(polyfit_modes,1);
% for i = 1:polyfit_modes
%     iter_list = 1:polyfit_modes;
%     iter_list = iter_list(iter_list ~= i);
%     L{i} = (high_res_low_eig(i)/low_eig(i)) * low_Minv_K;
%     for j = 1:(polyfit_modes-1)
%         L{i} = L{i} * (low_Minv_K - speye(size(low_res.M)) * low_eig(iter_list(j)))...
%             /(low_eig(i)-low_eig(iter_list(j)));
%     end
% end
% 
% Minv_K_new = L{1};
% for i = 2:polyfit_modes
%     Minv_K_new = Minv_K_new + L{i};
% end
hold on
lowdiag_new = lowV \ Minv_K_new * lowV;
lowdiag_new = diag(lowdiag_new);
% plot((lowdiag_new));
plot((lowdiag_new(1:plot_modes)));
plot((highD(1:plot_modes)));
legend('orignal low res', 'poly', 'high res')
% loweig_new = sort(-diag(lowdiag_new));
% [pV,pD] = eigs(Minv_K_new,12,'sm');
% [pV,pD] = eig(full(Minv_K_new));
% plowD = sort(diag(pD));
% plowD(1:polyfit_modes)

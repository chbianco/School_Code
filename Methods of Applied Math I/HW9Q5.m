%% Problem 5: SVD-Based Image Compression
% ENM 5200 HW9

clear; clc; close all;

%% Define the image matrix
A = [52 55 61 66 70 61 64 73 63 54;
     63 59 55 90 109 85 69 72 67 70;
     62 59 68 113 144 104 66 73 70 69;
     63 58 71 122 154 106 70 69 68 68;
     67 61 68 104 126 88 68 70 71 68;
     79 65 60 70 77 68 58 75 80 80;
     85 71 64 59 55 61 65 83 88 87;
     87 79 69 68 65 76 78 94 100 99;
     87 79 72 69 68 76 78 91 95 96;
     96 85 74 68 68 79 79 87 93 91];

%% Part 1: Compute the SVD
[U, S, V] = svd(A);
singular_values = diag(S);

fprintf('=== Part 1: Singular Values (descending order) ===\n');
for k = 1:length(singular_values)
    fprintf('  sigma_%d = %10.4f\n', k, singular_values(k));
end
fprintf('\n');

%% Part 2: Rank-1 Approximation
sigma1 = S(1,1);
u1 = U(:,1);
v1 = V(:,1);
A1 = sigma1 * (u1 * v1');

fprintf('=== Part 2: Rank-1 Approximation A1 ===\n');
disp(A1);

%% Part 3: Rank-2 Approximation
sigma2 = S(2,2);
u2 = U(:,2);
v2 = V(:,2);
A2 = A1 + sigma2 * (u2 * v2');

fprintf('=== Part 3: Rank-2 Approximation A2 ===\n');
disp(A2);

%% Part 4: Visualize and Quantify Errors

% Frobenius norm errors
err1_frob = norm(A - A1, 'fro');
err2_frob = norm(A - A2, 'fro');
A_frob    = norm(A, 'fro');

% Theoretical check: ||A - A_k||_F = sqrt(sum of remaining sigma_i^2)
err1_theory = sqrt(sum(singular_values(2:end).^2));
err2_theory = sqrt(sum(singular_values(3:end).^2));

% Spectral norm errors (largest discarded singular value)
err1_spec = norm(A - A1, 2);
err2_spec = norm(A - A2, 2);

% Energy captured (fraction of squared singular values retained)
energy1 = singular_values(1)^2 / sum(singular_values.^2);
energy2 = sum(singular_values(1:2).^2) / sum(singular_values.^2);

fprintf('=== Part 4: Error Analysis ===\n');
fprintf('||A||_F                = %8.4f\n', A_frob);
fprintf('||A - A1||_F (computed)= %8.4f\n', err1_frob);
fprintf('||A - A1||_F (theory)  = %8.4f  (= sqrt(sum_{i>=2} sigma_i^2))\n', err1_theory);
fprintf('||A - A2||_F (computed)= %8.4f\n', err2_frob);
fprintf('||A - A2||_F (theory)  = %8.4f  (= sqrt(sum_{i>=3} sigma_i^2))\n', err2_theory);
fprintf('||A - A1||_2           = %8.4f  (= sigma_2 = %.4f)\n', err1_spec, singular_values(2));
fprintf('||A - A2||_2           = %8.4f  (= sigma_3 = %.4f)\n', err2_spec, singular_values(3));
fprintf('Relative Frobenius error (rank-1): %.2f%%\n', 100*err1_frob/A_frob);
fprintf('Relative Frobenius error (rank-2): %.2f%%\n', 100*err2_frob/A_frob);
fprintf('Energy captured by rank-1: %.2f%%\n', 100*energy1);
fprintf('Energy captured by rank-2: %.2f%%\n', 100*energy2);

% Visualization: original, A1, A2, and error maps
figure('Position', [100 100 1200 700]);

subplot(2,3,1);
imagesc(A); colormap(gray); colorbar; axis image;
title('Original A'); caxis([min(A(:)) max(A(:))]);

subplot(2,3,2);
imagesc(A1); colormap(gray); colorbar; axis image;
title('Rank-1 Approximation A_1'); caxis([min(A(:)) max(A(:))]);

subplot(2,3,3);
imagesc(A2); colormap(gray); colorbar; axis image;
title('Rank-2 Approximation A_2'); caxis([min(A(:)) max(A(:))]);

subplot(2,3,4);
semilogy(singular_values, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Index k'); ylabel('$\sigma_k$ (log scale)');
title('Singular Value Spectrum'); grid on;

subplot(2,3,5);
imagesc(abs(A - A1)); colormap(gray); colorbar; axis image;
title('|A - A_1| (error map)');

subplot(2,3,6);
imagesc(abs(A - A2)); colormap(gray); colorbar; axis image;
title('|A - A_2| (error map)');

sgtitle('SVD-Based Image Compression');
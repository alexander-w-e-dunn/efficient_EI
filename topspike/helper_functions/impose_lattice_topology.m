function J_out = impose_lattice_topology(J_in, K)
% Enforce lattice topology and preserve input strength (weighted in-degree)

    N = size(J_in, 1);
    J_out = zeros(N);
    J_mask = zeros(N);  % binary mask of which connections are kept

    for i = 1:N
        for j = 1:K/2
            left  = mod(i - j - 1, N) + 1;
            right = mod(i + j - 1, N) + 1;

            J_mask(i, left)  = 1;
            J_mask(i, right) = 1;
        end
    end

    % Apply the mask
    J_out = J_in .* J_mask;

    % Renormalise column sums (preserve input to each neuron)
    in_strength_original = sum(J_in, 1);   % column sum = input strength
    in_strength_new = sum(J_out, 1);

    for j = 1:N
        if in_strength_new(j) > 0
            J_out(:,j) = J_out(:,j) * (in_strength_original(j) / in_strength_new(j));
        end
    end
end

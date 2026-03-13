function d_out = snellRefract(d_in, n_hat, n1, n2)
% Vector form of Snell's law.
% d_in    unit incident direction
% n_hat   unit surface normal pointing INTO the incident medium
% Returns d_out (unit refracted direction), or [] for TIR.

    d_in  = d_in  / norm(d_in);
    n_hat = n_hat / norm(n_hat);
    ratio = n1 / n2;
    cosI  = -dot(n_hat, d_in);
    sin2T = ratio^2 * (1 - cosI^2);

    if sin2T > 1.0   % total internal reflection
        d_out = [];
        return;
    end

    cosT  = sqrt(1 - sin2T);
    d_out = ratio * d_in + (ratio * cosI - cosT) * n_hat;
    d_out = d_out / norm(d_out);
end
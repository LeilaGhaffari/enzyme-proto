function dNdX = grad_basis_linear(xi)

    % derivatives of shape functions with respect to xi
    dN1_dxi = -0.125*(1-xi(2))*(1-xi(3));
    dN2_dxi = -dN1_dxi;
    dN3_dxi = 0.125*(1+xi(2))*(1-xi(3));
    dN4_dxi = -dN3_dxi;
    dN5_dxi = -0.125*(1-xi(2))*(1+xi(3));
    dN6_dxi = -dN5_dxi;
    dN7_dxi = 0.125*(1+xi(2))*(1+xi(3));
    dN8_dxi = -dN7_dxi;

    % derivatives of shape functions with respect to eta
    dN1_deta = -0.125*(1-xi(1))*(1-xi(3));
    dN2_deta = -0.125*(1+xi(1))*(1-xi(3));
    dN3_deta = -dN2_deta;
    dN4_deta = -dN1_deta;
    dN5_deta = -0.125*(1-xi(1))*(1+xi(3));
    dN6_deta = -0.125*(1+xi(1))*(1+xi(3));
    dN7_deta = -dN6_deta;
    dN8_deta = -dN5_deta;

    % derivatives of shape functions with respect to zeta
    dN1_dzeta = -0.125*(1-xi(1))*(1-xi(2));
    dN2_dzeta = -0.125*(1+xi(1))*(1-xi(2));
    dN3_dzeta = -0.125*(1+xi(1))*(1+xi(2));
    dN4_dzeta = -0.125*(1-xi(1))*(1+xi(2));
    dN5_dzeta = -dN1_dzeta;
    dN6_dzeta = -dN2_dzeta;
    dN7_dzeta = -dN3_dzeta;
    dN8_dzeta = -dN4_dzeta;

    % Populate dNdX
    dNdX = [dN1_dxi dN2_dxi dN3_dxi dN4_dxi dN5_dxi dN6_dxi dN7_dxi dN8_dxi;
            dN1_deta dN2_deta dN3_deta dN4_deta dN5_deta dN6_deta dN7_deta dN8_deta;
            dN1_dzeta dN2_dzeta dN3_dzeta dN4_dzeta dN5_dzeta dN6_dzeta dN7_dzeta dN8_dzeta];
end

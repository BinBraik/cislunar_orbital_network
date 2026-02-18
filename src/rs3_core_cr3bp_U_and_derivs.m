function pot = rs3_core_cr3bp_U_and_derivs(x,y,mu)
%RS3_CORE_CR3BP_U_AND_DERIVS  Potential U and derivatives for planar CR3BP.
% Reused from baseline nested cr3bp_U_and_derivs.
    x1 = x + mu;      y1 = y;
    x2 = x - 1 + mu;  y2 = y;
    r1_2 = x1.^2 + y1.^2;   r1 = sqrt(r1_2);
    r2_2 = x2.^2 + y2.^2;   r2 = sqrt(r2_2);
    r1_3 = r1_2 .* r1;
    r2_3 = r2_2 .* r2;
    r1_5 = r1_2 .* r1_3;           % r^5 = r^2 * r^3 (BUG FIX: was r^7)
    r2_5 = r2_2 .* r2_3;           % r^5 = r^2 * r^3 (BUG FIX: was r^7)

    pot.U  = 0.5*(x.^2 + y.^2) + (1-mu)./r1 + mu./r2;
    pot.Ux = x - (1-mu).*x1./r1_3 - mu.*x2./r2_3;
    pot.Uy = y - (1-mu).*y1./r1_3 - mu.*y2./r2_3;

    pot.Uxx = 1 + (1-mu).*(3*x1.^2 - r1_2)./r1_5 + mu.*(3*x2.^2 - r2_2)./r2_5;
    pot.Uyy = 1 + (1-mu).*(3*y1.^2 - r1_2)./r1_5 + mu.*(3*y2.^2 - r2_2)./r2_5;
    pot.Uxy =     (1-mu).*(3*x1.*y1)./r1_5     + mu.*(3*x2.*y2)./r2_5;
end

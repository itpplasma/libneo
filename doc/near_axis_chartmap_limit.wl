(* ::Package:: *)

(* Analytic near-axis limit of the Boozer chartmap (rho = Sqrt[s]).
   Derives the position, covariant basis (chartmap metric), Jacobian and the
   pseudo-Cartesian regularization of the flux->Cartesian map as rho->0, to
   implement the exact s=0 slice of the chartmap export without evaluating VMEC
   surface harmonics at the axis (no surface exists there).

   Stellarator-symmetric VMEC map with near-axis (VMEC / Garren-Boozer)
   regularity R_m(rho) ~ rho^|m|: only m=0 survives at the axis (the axis curve),
   m=1 gives the limiting ellipse. Position in cylindrical (R, phi, Z), then
   Cartesian x = R Cos[phi], y = R Sin[phi], z = Z. *)

ClearAll["Global`*"];

(* Near-axis scalar expansions to O(rho^2). Coefficients are functions of zeta
   ze (the toroidal coordinate). phi is the cylindrical toroidal angle, regular
   at the axis (phi -> p0[ze]); only the poloidal structure carries rho^|m|. *)
Rr = r0[ze] + rho (r1c[ze] Cos[th] + r1s[ze] Sin[th])
   + rho^2 (r20[ze] + r2c[ze] Cos[2 th] + r2s[ze] Sin[2 th]);
Zz = z0[ze] + rho (z1c[ze] Cos[th] + z1s[ze] Sin[th])
   + rho^2 (z20[ze] + z2c[ze] Cos[2 th] + z2s[ze] Sin[2 th]);
ph = p0[ze] + rho (p1c[ze] Cos[th] + p1s[ze] Sin[th]);

xvec = {Rr Cos[ph], Rr Sin[ph], Zz};

Print["=== Position x = (R Cos phi, R Sin phi, Z), near-axis expansion ==="];
Print["x  = ", Series[xvec[[1]], {rho, 0, 1}]];

(* Covariant basis columns dx/d(rho,theta,zeta) *)
eRho = D[xvec, rho];
eThe = D[xvec, th];
eZet = D[xvec, ze];

aRho = Simplify[eRho /. rho -> 0];
aThe = Simplify[eThe /. rho -> 0];
aZet = Simplify[eZet /. rho -> 0];
Print["\n=== Covariant basis at the axis (rho=0) ==="];
Print["dx/drho |0   = ", aRho, "   (m=1 ellipse vectors; NONZERO)"];
Print["dx/dtheta|0  = ", aThe, "   (vanishes ~rho: degenerate)"];
Print["dx/dzeta |0  = ", aZet, "   (axis tangent)"];

(* Metric and Jacobian leading behaviour *)
gTT = Simplify[eThe . eThe];
gRR = Simplify[eRho . eRho];
Jpolar = Det[Transpose[{eRho, eThe, eZet}]];
Print["\n=== Metric / Jacobian leading order in rho ==="];
Print["g_rho_rho(0) = ", Simplify[gRR /. rho -> 0], "   (finite)"];
Print["g_theta_theta = ", Normal[Series[gTT, {rho, 0, 2}]], "   (O(rho^2) -> 0)"];
Print["det[J_polar]  = ", Simplify[Normal[Series[Jpolar, {rho, 0, 1}]]],
   "   (O(rho): sqrt(g) ~ rho, vanishes at axis)"];

(* Pseudo-Cartesian chart X = rho Cos[th], Y = rho Sin[th] : regular at axis.
   Rewrite the O(rho^2) position with rho Cos[th]=Xv, rho Sin[th]=Yv,
   rho^2 = Xv^2+Yv^2, rho^2 Cos[2th] = Xv^2-Yv^2, rho^2 Sin[2th] = 2 Xv Yv. *)
sub = {rho Cos[th] -> Xv, rho Sin[th] -> Yv, rho^2 -> Xv^2 + Yv^2,
   rho^2 Cos[2 th] -> Xv^2 - Yv^2, rho^2 Sin[2 th] -> 2 Xv Yv};
RrC = r0[ze] + (r1c[ze] Xv + r1s[ze] Yv)
   + r20[ze] (Xv^2 + Yv^2) + r2c[ze] (Xv^2 - Yv^2) + r2s[ze] (2 Xv Yv);
ZzC = z0[ze] + (z1c[ze] Xv + z1s[ze] Yv)
   + z20[ze] (Xv^2 + Yv^2) + z2c[ze] (Xv^2 - Yv^2) + z2s[ze] (2 Xv Yv);
phC = p0[ze] + (p1c[ze] Xv + p1s[ze] Yv);
xC = {RrC Cos[phC], RrC Sin[phC], ZzC};
eX = D[xC, Xv] /. {Xv -> 0, Yv -> 0};
eY = D[xC, Yv] /. {Xv -> 0, Yv -> 0};
eZc = D[xC, ze] /. {Xv -> 0, Yv -> 0};
Print["\n=== Pseudo-Cartesian (X,Y) basis at axis (regular) ==="];
Print["dx/dX|0 = ", Simplify[eX]];
Print["dx/dY|0 = ", Simplify[eY]];
Print["det[J_(X,Y,ze)]|0 = ", Simplify[Det[Transpose[{eX, eY, eZc}]]],
   "   (finite, nonzero -> chart regular through the axis)"];

Print["\n=== Identification with VMEC data ==="];
Print["axis (m=0): R0(ze)=p? no -> r0 = Sum_n raxis_cc(n) Cos[n nfp ze],"];
Print["            Z0(ze) = Sum_n zaxis_cs(n) Sin[n nfp ze]  (exact, from wout)."];
Print["x_axis = (r0 Cos[p0], r0 Sin[p0], z0),  p0 = cyl toroidal angle limit."];
Print["ellipse (dx/drho|0): r1c,r1s,z1c,z1s = the m=1 harmonics' rho-slope."];
Print["A_phi(0)=0 (no enclosed flux); |B|,B_theta,B_phi from healed m=0 at rho=0."];
Print["\nDONE"];

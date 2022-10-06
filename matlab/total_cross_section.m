% \brief Calculate total cross section for given parameters.
%
% Calculate the total cross section of a particle with energy E
% colliding with a particle at rest, using the model equation
%
% \sigma = (A_5 + A_2 / ((A_4 - A_3*E)^2 + 1)) / (E*(\exp(A_1 / \sqrt(E)) - 1))
%
% Cross sections for different energies can be calculated with one call,
% while multiple parameter sets require multiple calls.
%
% Example:
%
% E = 1:200;
% cross_section = total_cross_section(E, [45.95, 50200, 1.368e-2, 1.076, 409]);
%
% to calculate a range of cross sections. Using
%
% cross_section = total_cross_section(E, 'dt');
%
% should give the same result.
%
% input:
% ------
% E: floating point number, the energy in keV of the particle.
%   Might also be an array.
% A: vector of five numbers, the parameters to use for the model
%   or a string, determining one of the predefined values (not case
%   sensitive).
%   DDa: deuterium - deuterium channel a (T+p)
%   DDb: deuterium - deuterium channel b (He3 + n)
%   DHe or DHe3: deuterium - helium-3
%   DT: deuterium - tritium
%   THe or THe3: tritium - helium-3, all three chanels
%   TT: tritium - tritium
%   If a char/string is given for A, and no match is found, exectuion
%   of the function is stopped with an error message.
%
% output:
% -------
% cross_section: cross section according to the model for the given
%   energy. Will have the same shape as input E.
function cross_section = total_cross_section(E, A)
  if ischar(A)
    A = lower(A);
    switch(A)
    case 'dda'
      A = [46.097, 372, 4.36e-4, 1.220, 0.0];
    case 'ddb'
      A = [47.88, 482, 3.08e-4, 1.177, 0.0];
    case {'dhe', 'dhe3'}
      A = [89.27, 25900, 3.98e-3, 1.297, 647];
    case 'dt'
      A = [45.95, 50200, 1.368e-2, 1.076, 409];
    case {'the', 'the3'}
      A = [123.1, 11250, 0.0, 0.0, 0.0];
    case 'tt'
      A = [38.39, 448, 1.02e-3, 2.09, 0.0];
    otherwise
      error(['Unknown string value for parameter A = ', A]);
    end
  end

  cross_section = (A(5) + A(2) ./ ((A(4) - A(3)*E).^2 + 1)) ./ (E.*(exp(A(1) ./ sqrt(E)) - 1));
end

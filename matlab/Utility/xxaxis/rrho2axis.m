function rrho2axis(a, b, r, rho2, delta)
%##########################################################################
%function rrho2axis(a, b, r, rho2, delta)
%##########################################################################
% description:
%--------------------------------------------------------------------------
% converts the 2nd x-axis where the 1st axis (a) is given on r to rho. The
% variables r and rho give the corresponding table to calculate between the
% axis where linear interpolation and extrapolation is used. With delta one
% can change the uniform rho spacing.
%##########################################################################
% input:
%--------------------------------------------------------------------------
% a     ... 1st axis (on r)
% b     ... 2nd axis (will be converted to rho)
% r     ... r values
% rho2  ... rho values
% delta ... uniform rho spacing
%##########################################################################

%author:   Philipp Ulbl
%created:  04.03.2020

    %get ticks of r axis
    ticka_r = get(a, 'XTick');
    %convert ticks to rho
    ticka_rho = interp1(r, rho2, ticka_r, 'linear', 'extrap');
    
    %create ticks for rho axis by rounding converted r axis
    tickb_rho_raw = round(ticka_rho, 2);
    %tickb_rho = tickb_rho(tickb_rho<=1); %maybe exclude rho > 1 ?
    
    %get spacing of uniform rho ticks by parameter or calculate
    if(nargin < 5 || isempty(delta))
        delta_psi = ceil((ticka_rho(end)-ticka_rho(end-1)) * 100) / 100;
    else
        delta_psi = delta;
    end
    
    %create spacing in rho with 0.95 and 1 included
    tickb_rho = tickb_rho_raw(1):delta_psi:min(0.92, tickb_rho_raw(end));
    if(tickb_rho_raw(end) >= 0.95)
        if(delta_psi <= 0.05)
            tickb_rho = [tickb_rho, 0.95];
        end
        tickb_rho = [tickb_rho, 1];
        if(delta_psi <= 0.05)
            tickb_rho = [tickb_rho, 1.03:delta_psi:tickb_rho_raw(end)];
        end
    end
    
    %recalculate r ticks for rho axis
    tickb_r = interp1(rho2, r, tickb_rho, 'linear', 'extrap');

    %set r ticks for rho axis
    set(b,'XTick', tickb_r)
    %and label them by rho
    set(b,'XTickLabel', arrayfun(@num2str, tickb_rho, 'UniformOutput', false));
    %set limit to the same as r axis
    xlim(get(a, 'Xlim'))
    %label as rho poloidal
    xlabel('\rho_{pol}^2')
end

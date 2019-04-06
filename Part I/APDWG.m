function [ neff_TE , neff_TM , Ey_TE , Ex_TM, Hy_TM ] = APDWG( wl , h , n1 , n2 , n3 , x )

% ADPWG -- Asymmetric Planar Dielectric WaveGuide (Original version Mar 2013)
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
%                        March 2018 version
%  Returns the Ex component instead of Hy component for TM modes      
%  It has been improved (lines 190, 196) to avoid the occurence of Infs
%  These are the only differences compared to the original 2013 version 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% ==================================================================
% | Inputs:
% ==================================================================
%   wl : wavelength, e.g. wl=1.55;
%   h  : thickness of guiding-layer (slab), e.g. h=0.22;
%   n1 : guiding-layer refractive index, e.g. n1=3.45;
%   n2 : substrate refractive index (n2<n1), e.g. n2=1.45;
%   n3 : cladding refractive index (n3<n1), e.g. n1=1.00;
%   x  : width vector over which to calculate values of Ey and Hy 
%        for TE/TM, respectively. E.g. x=linspace(-2*h,+3*h,1000);
% ==================================================================
% | Outputs:
% ==================================================================
%   neff_TE : effective indices of all guided TE-modes
%   neff_TM : effective indices of all guided TM-modes
%   Ey_TE   : profile(s) of Ey component, over the x-vector, 
%             for all guided TE-mode(s)
%   Ex_TM   : profile(s) of Ex component, over the x-vector, 
%             for all guided TM-mode(s)
% ------------------------------------------------------------------
% | Notes:
% ------------------------------------------------------------------
% (1) Inputs wl, h and x should have the same units, e.g. um
% (2) If x is omitted (or x=[];) the Ey_TE & Hy_TM profiles are not 
%     calculated.
% (3) If n3==n2, then we have an SPDWG. Please remember that, unlike
%     the SPDWG, the APDWG has a non-zero cutoff frequency (or width)
%     even for the fundamental mode (TE & TM).
% (4) For very large V-numbers (many, many modes), or for modes too
%     close to their cut-off (h or wl), the root-solving algorithm 
%     employed by APDWG might miss some roots of the Characteristic
%     Equation. In such an event, one solution would be to increase 
%     the n-eff search-range discretization from its default value,
%     Nne=50 (search the code below). Anyway, such situations should
%     be avoided, so its better to work at lower V-numbers.
% (5) When called with no output arguments, this APDWG function 
%     plots the Characteristic Equation as a function of n_eff, and 
%     the mode-profiles, for TE & TM (x must be supplied).
% 
% Alexandros Pitilakis, Thessaloniki, March 2013

%Test Input Vars:
if nargin == 0
    n1 = 3.45; % Index : Guiding Layer (Silicon)
    n2 = 1.45; % Index : Substratre (Oxide)
    n3 = 1.00; % Index : Cladding (Air)
    wl = 1.55; %[um] Wavelength
    h  = 0.22; %[um] Slab Thickness
    x  = linspace( -2*h , 2*h , 1000 ); %[um] transverse-coordinate vector
elseif nargin == 5
    x = []; % This means that no Hy,Ey output field-profiles were requested
elseif nargin == 4
    n3 = n2; % SDPWG
    x = [];
end

% Error-Checks:
if any( ~isreal( [n1 n2 n3 wl h] ) )
    error( ' ## APDWG: Inputs are strictly REAL!' ); 
end
if n2>=n1 || n3>=n1
    error( ' ## APDWG: n1>n2 AND n1>n3! (1=Guiding, 2=Substrate, 3=Cladding)' ); 
end
if ~isempty(x) && ( min(x)>0 || max(x)<h )
    error( ' ## APDWG: x-vector must be equal-to or wider-than the [0,h] x-region!' )
end
if isempty(x) && nargout == 0
    error( ' ## APDWG: x-vector must be supplied for no-output plots!' )
end

% Simulation Info
if nargout == 0
    disp(    '  -------> Launching APDWG <------- ' );
    fprintf( '  ** Wavelength     (wl) = %6.4fum\n' , wl );
    fprintf( '  ** Slab Thickness  (h) = %6.4fum\n' , h  );
    fprintf( '  ** Guiding Layer  (n1) = %6.4f\n' , n1  );
    fprintf( '  ** Substrate      (n2) = %6.4f\n' , n2  );
    fprintf( '  ** Cladding       (n3) = %6.4f\n' , n3  );
    disp(    '  --------------------------------- ' );
end

%-------------------------------------------------------------------------
% Mode-Number Analysis
%-------------------------------------------------------------------------

k0 = 2*pi / wl; % [um^-1] Wavenumber, free-space

xi=(n2^2-n3^2)/(n1^2-n2^2); %APDWG "asymmetry" factor (0==SPDWG)
V = k0*h/2*sqrt( n1^2-n2^2 ); %V-number, same definition as for SPDWG

VasosTE = atan( sqrt(xi) )/2; % V-number "asymmetry-offset" for TE-modes
VasosTM = atan( (n1/n3)^2 * sqrt(xi) )/2; % V-number "asymmetry-offset" for TM-modes

nModeTE = ceil( (V-VasosTE)/(pi/2) ); %Number of TE-modes
nModeTM = ceil( (V-VasosTM)/(pi/2) ); %Number of TM-modes
clear VasosTE VasosTM xi V

%-------------------------------------------------------------------------
% Infinites/Singularities/Discontinuities of the Characteristic Equation
%-------------------------------------------------------------------------

% From close examination of the particular transcedental Char.Eqs of the
% APDWG, it can be observed that all the n_eff roots can be found between
% consecutive discontinuities, or other "well defined" values of refr indx.
% The discontinuitis or infinities occur from "tan(x)" and "1/x" terms 
% of the Ch.Eq. when x=pi/2 and 0, respectively.

e1 = n1^2;
e2 = n2^2;
e3 = n3^2;

%For both TE/TM : Terms from tan(kh)=Inf
ms = 1:2:floor( sqrt(e1-e2)*h*k0*2/pi ); % odd-multiples
nInfs_tan = sqrt( e1 - (ms*pi/2/h/k0).^2 );

%For TE: Terms from k^2=g*d
nInfs1 = sqrt( (-e1^2+e2*e3)/(-2*e1+e2+e3) );

%For TM: Terms from k^2=g'*d'
A = e2*e3/e1^2; 
nInfs2 = sqrt( roots( [ A^2-1 , -2*e1*A^2+e2+e3 , e1^2*A^2-e2*e3 ] ) );
nInfs2 = nInfs2( nInfs2>n2 & nInfs2<n1 )';

%Altogether now:
nInfsTE = sort( [ nInfs1(:) ; nInfs_tan(:) ] )';
nInfsTM = sort( [ nInfs2(:) ; nInfs_tan(:) ] )';

clear e1 e2 e3 A ms nInfs1 nInfs2 nInfs_tan

% ========================================================================
% Solve Characteristic Equation for n_effective
% ========================================================================
for kT = 1 : 2 % {1,2}={TE,TM} modes

    % This is the number of samples in the specified neff-range where 
    % we'll look for roots using the interp1 function. Read on.
    Nne = 50;
        
    % Initialize params, depending on TE,TM
    switch kT
        case 1, 
            N = nModeTE; % Number of modes expected
            nI = nInfsTE;  % The n-effs of the "Inf" values
            neff_TE = NaN*ones(1,N); % Holds the n-effs of the TE-modes
            neTE = NaN*ones(N,Nne); % Aux: Holds the n-effs for ChEq plots
            XETE = NaN*ones(N,Nne); % Aux: Values of ChEq for the neTE above
        case 2, 
            N = nModeTM; 
            nI = nInfsTM; 
            neff_TM = NaN*ones(1,N);
            neTM = NaN*ones(N,Nne);
            XETM = NaN*ones(N,Nne);
    end

    % Check the number of pairs of infinities. If they're less than the
    % number of modes expected, add the substrate value as the last limit value.
    if length(nI)-1 < N
        nI = sort( [ max( [n2,n3,1] ) , nI ] );
    end
    
    % Scan modes
    for kn = 1 : N

        % Define a discrete n_eff range, between the limits given above.
        % Mutliply limit values with 1+/-eps (the machine precision) 
        % to get rid of numerical errors.
        ne = linspace( nI(kn)*(1+eps) , nI(kn+1)*(1-eps)  , Nne ); 

        %Characteristic Equation Parameters
        k  = k0 * sqrt( n1^2 - ne.^2 ); %common for TE/TM
        
        g  = k0 * sqrt( ne.^2 - n2^2 ); %TE modes
        d  = k0 * sqrt( ne.^2 - n3^2 ); %TE modes
        
        gt = g * (n1/n2)^2;             %TM modes
        dt = d * (n1/n3)^2;             %TM mode

        %Form the Chacteristic Eqs & Solve Numerically w Cubic Interpolation
        switch kT
            case 1, 
                XE  = tan( k * h ) - k.*(g +d )./(k.^2-g.*d  ); %TE
                XE(isinf(XE)) = sign(XE(isinf(XE)))/eps; % Replace Infs with something really large, but finite. Infs annoy INTERP1 function.             
                neff_TE(kn) = interp1( XE , ne , 0 , 'pchip');
                neTE(kn,:) = ne;
                XETE(kn,:) = XE;
            case 2, 
                XE = tan( k * h ) - k.*(gt+dt)./(k.^2-gt.*dt); %TM
                XE(isinf(XE)) = sign(XE(isinf(XE)))/eps; % Replace Infs with something really large, but finite. Infs annoy INTERP1 function.
                neff_TM(kn) = interp1( XE , ne , 0 , 'pchip' );
                neTM(kn,:) = ne;
                XETM(kn,:) = XE;
        end

    end; clear kn 

end; clear kT ne XE

Nmodes_TE = length( neff_TE );
Nmodes_TM = length( neff_TM );

%Cross-Check Results for Number of modes / Analytical vs Numerical
if length( neff_TE ) ~= nModeTE
    error( ' ## APDWG: Analytical vs Numerical Number of TE modes do not match! Increse Nne (?)' );
elseif length( neff_TM ) ~= nModeTM
    error( ' ## APDWG: Analytical vs Numerical Number of TM modes do not match! Increse Nne (?)' );
end

% If some modes are not found...
if Nmodes_TE == 0 
    disp(    ' ## APDWG: No TE modes supported! Increase V-number (?).' );
    fprintf( '    [wl_um, h_um, n1, n2, n3] = [%6.4f, %6.4f, %4.2f, %4.2f, %4.2f]\n' , wl,h,n1,n2,n3 );
end
if Nmodes_TM == 0 
    disp( ' ## APDWG: No TM modes supported! Increase V-number (?).' );
    fprintf( '    [wl_um, h_um, n1, n2, n3] = [%6.4f, %6.4f, %4.2f, %4.2f, %4.2f]\n' , wl,h,n1,n2,n3 );
end
if Nmodes_TE == 0 && Nmodes_TM == 0 , return; end


% ========================================================================
% Calc Field Distributions
% ========================================================================

% If no output-fields were requested...
if isempty( x )
    if nargout > 2 
        Ey_TE = [];
        Ex_TM = [];
        Hy_TM = [];
    end
    
    %The ChEq roots (neffs of supported modes) are sorted in ascending
    %order. We, usually, prefer to call #1 the fundamental mode, who in the
    %above ordering comes last, hence the flipping/reversion
    neff_TE = fliplr( neff_TE );
    neff_TM = fliplr( neff_TM );    
    
    return
end

% Translate x-vector, to have its center at the middle of the guide
%x = x + h/2;       % commented out by EEK

% Region-Separation boolean x-vectors for Ey(x) and Hy(x).
is1 = x>=0 & x<=h ; %Guiding Layer
is2 = x<0         ; %Substrate
is3 = x>h         ; %Cladding

% TE-Modes
Ey_TE = NaN * zeros( Nmodes_TE , length(x) );
for kkm = 1 : Nmodes_TE % Scan modes
    neff = neff_TE(kkm);
    k  = k0 * sqrt( n1^2 - neff^2 );
    g  = k0 * sqrt( neff^2 - n2^2 );
    d  = k0 * sqrt( neff^2 - n3^2 );
    phi = atan( g / k );
    A1 = 1 ;
    A2 = A1 * cos( phi );
    A3 = A1 * cos( k*h -phi );
    Ey_TE( kkm , is1 ) = A1 * cos( k * x(is1) - phi );
    Ey_TE( kkm , is2 ) = A2 * exp( g * x(is2) );
    Ey_TE( kkm , is3 ) = A3 * exp( -d *(x(is3)-h) );
end

%The Refractive Indices across the x-vector cross-section
% // These are used for the TM-modes' Ex(x)-calculation
ns = NaN*zeros(1,length(x));
ns( is1 ) = n1;
ns( is2 ) = n2;
ns( is3 ) = n3;

% TM-Modes
Hy = NaN * zeros( 1 , length(x) );
Hy_TM = NaN * zeros( Nmodes_TM , length(x) );
Ex_TM = NaN * zeros( Nmodes_TM , length(x) );
for kkm = 1 : Nmodes_TM % Scan modes
    neff = neff_TM(kkm);
    k  = k0 * sqrt( n1^2 - neff^2 );
    g  = k0 * sqrt( neff^2 - n2^2 );
    d  = k0 * sqrt( neff^2 - n3^2 );
    phi = atan( g / k  * (n1/n2)^2 );
    A1 = 1 ;
    A2 = A1 * cos( phi );
    A3 = A1 * cos( k*h -phi );
    Hy( is1 ) = A1 * cos( k * x(is1) - phi );
    Hy( is2 ) = A2 * exp( g * x(is2) );
    Hy( is3 ) = A3 * exp( -d *(x(is3)-h) );
    Hy_TM( kkm , : ) = Hy ./ max(abs(Hy)) ;
    Ex = Hy ./ ns.^2 ; %Ex is proportional (not equal) to Hy/n^2, from TM-wave Maxwell's eqs.
    Ex_TM( kkm , : ) = Ex ./ max(abs(Ex)) ;
  
end

% ========================================================================
% Flip Result Order!
% ========================================================================

%The roots (neffs of supported modes) and fields are sorted in ascending
%order. We, usually, prefer to call #1 the fundamental mode, who in the
%above ordering comes last, hence the flipping/reversion
neff_TE = fliplr( neff_TE );
neff_TM = fliplr( neff_TM );
Ey_TE = flipud( Ey_TE );
Hy_TM = flipud( Hy_TM );
Ex_TM = flipud( Ex_TM );

%-------------------------------------------------------------------------
% Plots
%-------------------------------------------------------------------------
if nargout == 0

    if Nmodes_TE > 10 || Nmodes_TM > 10
        fprintf( ' #### Too many modes to plot ( %dxTE & %dxTM )... You do it...\n' , ...
            Nmodes_TE , Nmodes_TM );
        return
    end

    % Prepare Figure
    FigLbl1 = 'APDWG' ;
    if isempty( findobj( 'Name',FigLbl1) )
        figure; set(gcf,'NumberTitle','off','Name',FigLbl1); hold on;
    else
        figure( findobj( 'Name',FigLbl1) ); clf;
    end; clear FigLbl1

    % ....................................................................
    % Characteristic Equation Solving -> n_eff // TE
    % ....................................................................
    if Nmodes_TE > 0
        subplot(2,2,1); hold on
        ne = reshape( neTE' , Nmodes_TE*Nne , 1 );
        XE = reshape( XETE' , Nmodes_TE*Nne , 1 );
        plot( real(ne) , real(XE) , 'bs' ); hold on
        plot( real(ne) , zeros(size(ne)) , 'k' , 'LineWidth' , 2 );
        axis( [min(ne) max(ne) 10*[-1 1]] )
        axis on; grid on;
        for ku = 1:length(nInfsTE) % Plot Infinities of ChEq
            plot( nInfsTE(ku)*[1 1] , get(gca,'YLim') , 'b-' , 'LineWidth' , 2  )
        end; clear ku
        %Annotations
        str1 = 'function( n_{eff} ) = tan(\kappah) - \kappa(\gamma+\delta)/( \kappa^2 - \gamma\delta)';
        ylabel( str1 ); xlabel( 'n_{effective} = \beta / k_0' );
        title( sprintf( 'APDWG Char. Equation for TE-modes' ) );
        text( neff_TE(1) , -0.25 , sprintf( '\\uparrow neff(TE) ~ %6.4f' , neff_TE(1) ) , ...
            'FontWeight' , 'Bold' , 'FontSize' , 12 , 'VerticalAlignment' , 'Top' , ...
            'Color' , 'b' , 'EdgeColor' , 'b' , 'BackGroundColor' , 'w' )
    end

    % ....................................................................
    % Characteristic Equation Solving -> n_eff // TM
    % ....................................................................
    if Nmodes_TM > 0
        subplot(2,2,2); hold on
        ne = reshape( neTM' , Nmodes_TM*Nne , 1 );
        XE = reshape( XETM' , Nmodes_TM*Nne , 1 );
        plot( real(ne) , real(XE) , 'ro' ); hold on
        plot( real(ne) , zeros(size(ne)) , 'k' , 'LineWidth' , 2 );
        axis( [min(ne) max(ne) 10*[-1 1]] )
        axis on; grid on;
        for ku = 1:length(nInfsTM) % Plot Infinities of ChEq
            plot( nInfsTM(ku)*[1 1] , get(gca,'YLim') , 'r-' , 'LineWidth' , 2  )
        end; clear ku
        %Annotations
        str1 = 'function( n_{eff} ) = tan(\kappah) - \kappa(\gamma''+\delta'')/( \kappa^2 - \gamma''\delta'')';
        ylabel( str1 ); xlabel( 'n_{effective} = \beta / k_0' );
        title( sprintf( 'APDWG Char. Equation for TM-modes' ) );
        text( neff_TM(1) , +0.25 , sprintf( '\\downarrow neff(TM) ~ %6.4f' , neff_TM(1) ) , ...
            'FontWeight' , 'Bold' , 'FontSize' , 12 , 'VerticalAlignment' , 'Bottom' , ...
            'Color' , 'r' , 'EdgeColor' , 'r' , 'BackGroundColor' , 'w' )
    end

    % ....................................................................
    % Field-Profiles // TE
    % ....................................................................
    if Nmodes_TE > 0
        subplot(2,2,3); hold on; grid on;

        colsTE = winter( Nmodes_TE );
        Lws = fliplr( linspace( 1 , 3 , Nmodes_TE ) );
        Legs=cell(1,Nmodes_TE);
        for kmte = 1 : Nmodes_TE
            plot( x ,Ey_TE(kmte,:) , 'LineWidth' , Lws(kmte) , 'Color' , colsTE(kmte,:) , 'LineStyle' , '-' );
            Legs{kmte}=sprintf( 'TE%d' , kmte-1 );
        end
        legend( Legs ); clear Legs

        plot( x , 0*x , 'k' );
        if n3 == 1, cCla = 1; else cCla = 0.80; end
        if n2 ~= n3, cSub = 0.75; else cSub = cCla; end
        fill( [min(x) min(x) 0 0] , 1*[-1 1 1 -1] , cSub*[1 1 1] , 'FaceAlpha' , 0.5 ); %Substrate
        fill( [0 0 h h]           , 1*[-1 1 1 -1] , 0.50*[1 1 1] , 'FaceAlpha' , 0.5 ); %Guiding
        fill( [h h max(x) max(x)] , 1*[-1 1 1 -1] , cCla*[1 1 1] , 'FaceAlpha' , 0.5 ); %Cladding
        set( gca,'XLim',[min(x) max(x)], 'YLim' , [-A1 A1] );
        xlabel( 'x-axis [\mum]' ); ylabel( 'Amplitude E_y(x)' );
        title( 'TE mode(s)' )
    end

    % ....................................................................
    % Field-Profiles // TM  -- Older version (in comments) plots Ex_TM piecewise to display the discontinuity
    % ....................................................................
    if Nmodes_TM > 0
        subplot(2,2,4); hold on; grid on;

        colsTM = autumn(Nmodes_TM) ;
        Lws = fliplr( linspace( 1 , 3 , Nmodes_TM ) );
        Legs=cell(1,Nmodes_TM); 
        hp1=NaN*ones(1,Nmodes_TM);
        for kmtm = 1 : Nmodes_TM
            ib(1) = find( x<0 , 1 , 'last' );
            ib(2) = find( x<h , 1 , 'last' );
            hp1(kmtm)=plot( x((      1):ib(1)) , (Ex_TM(kmtm,(      1):ib(1))) , 'LineWidth' , Lws(kmtm) , ...
                'Color' , colsTM(kmtm,:) , 'LineStyle' , '-' );
            plot( x((ib(1)+1):ib(2)) , (Ex_TM(kmtm,(ib(1)+1):ib(2))) , 'LineWidth' , Lws(kmtm) , ...
                'Color' , colsTM(kmtm,:) , 'LineStyle' , '-' );
            plot( x((ib(2)+1):end  ) , (Ex_TM(kmtm,(ib(2)+1):end  )) , 'LineWidth' , Lws(kmtm) , ...
                'Color' , colsTM(kmtm,:) , 'LineStyle' , '-' );
            %plot( x ,Hy_TM(kmtm,:) , 'LineWidth' , Lws(kmtm) , 'Color' , colsTM(kmtm,:) , 'LineStyle' , '-' );
            Legs{kmtm}=sprintf( 'TM%d' , kmtm-1 );
        end
        %legend( hp1 , Legs ); clear Legs
        legend( Legs ); clear Legs

        plot( x , 0*x , 'k' );
        if n3 == 1, cCla = 1; else cCla = 0.80; end
        if n2 ~= n3, cSub = 0.75; else cSub = cCla; end
        fill( [min(x) min(x) 0 0] , 1*[-1 1 1 -1] , cSub*[1 1 1] , 'FaceAlpha' , 0.5 ); %Substrate
        fill( [0 0 h h]           , 1*[-1 1 1 -1] , 0.50*[1 1 1] , 'FaceAlpha' , 0.5 ); %Guiding
        fill( [h h max(x) max(x)] , 1*[-1 1 1 -1] , cCla*[1 1 1] , 'FaceAlpha' , 0.5 ); %Cladding
        set( gca,'XLim',[min(x) max(x)], 'YLim' , [-A1 A1] );
        xlabel( 'x-axis [\mum]' ); 
        ylabel( 'Amplitude E_x(x)' );
        %ylabel( 'Amplitude H_y(x)' );
        title( 'TM mode(s)')
    end

end








% evalPSDBBDD() - Subroutine to comupute FFT spectrum and obtain
% Power Spectral Density (PSD), both unfiltered (broadband, BB)
% and filtered in decidecade bands (DD sound pressure level).
% Called by 'TUNE' main.
% Use 'periodogram' with given FFT no. within a single snapshot, as
% many times as needed to fill each snapshot with given overlap.
% Use the given calibration data converted to Scale Factor in Pa.
% Output normalization to window total power is already accounted for 
% in function 'periodogram()', while frequency bin width normalization
% is added explicitely in DD computation.
% Input arguments: snapshot, FFT window, calib.SF (BB and DD), 
% DD frequencies, parameters. DD arguments are void if not required.
% Output are PSD matrixes in dB re 1 microPa^2/Hz with hi/avg/low values.
% Averaging on snapshots may be either arithmetic mean or median.
% Previously used zero padding for last FFT subinterval is no longer used:
% instead, last subinterval is slid back to fit in last snapshot portion.
% DD matrix provides band limits and indexes of f.vector.
% Draw BB or DD plot only if corresponding option is enabled in main.
% Debug mode requires key strokes between plots.
%
% Reference: ISO 18405:2017 3.1.3.13, 3.1.4.2; Lurton p.132 for Knudsen
%
% REVISION HISTORY
% Note: merged from earlier 'evalPSD', 'evalTOB' development subroutines.
% 3.0 - 27-mar-2024 - Merged from 'evalPSD26' and 'evalTOB27' subroutines.
%       Changed frequency intervals from one-third octave to decidecade.
%       BB elements are summed within DD band and divided by bandwidth,
%       then normalized to f bin width using srate/nfft factor.
%       Earlier subroutines used band averaging without normalization.
%       Added WAV figure handle and boolean proc.type as input param.
%       Input arguments are void if processing is not required.
%       Two output arguments - BB and DD matrixes - instead of one.
%       Requires main calling function 'UWNoiPro506main' or later.
% 3.1 - 26-nov-2024 - Fixed error due to return value not assigned when
%	not using BB option. Improved zeropad. Requires main v. 5.10.
% 3.2 - 27-nov-2024 - Alternate handling of last FFT subinterval:
%       instead of zeropad, slide back & use last full subinterval. Fixes
%       issues found at very low frequency (below about 10 Hz) as
%       interpolating zeropad artificially increases low-f components.
%       Last FFT subinterval now always overlaps with previous one.
%       Requires main v. 5.11.
% 3.3 - Removed zero padded subinterval normalization, as sliding back
%       maintains same subinterval length. DD PSD redefined as DD SPL,
%       as well as related setup parameters. Requires 'TUNE513main'.
% 3.4 - Aligned version with 'TUNE514main', minor code changes.
% 4.0 - Decidecade SPL is now a sum of components and not an average
%   	(removed division by bandwidth). Also provide full band power as 
%       sum over entire spectrum: was done before in main as wave rms. 
%       It allows weighting by SF frequency components instead of using 
%       mean SF, which was less accurate. Requires 'TUNE601main'.
%
% COPYRIGHT
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    any later version.
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See GNU General
%    Public License for more details  <https://www.gnu.org/licenses/>.
%
% AUTHORS
%   Silvano Buogo, Junio Fabrizio Borsani, Valentina Caradonna
%   <silvano.buogo@cnr.it>
%
% VERSION 4.0 | 02-apr-2025 

function [bbpxxavg,ddpxxavg,fbpxxavg] = ...
                evalPSDBBDD40(Snap,Win,BBSF,DDSF,DDmat,Param)

% Parameters:
% - Snap.Y   = snapshot array of samples
% - Snap.T   = snapshot time array in s
% - Win      = time window
% - SF       = calibration (scale factor in Pa) for BB frequencies
% - DDSF     = as above, for DD center frequencies
% - DDmat    = decidecade band matrix with center and edge freq.s
% - Param.sr = sample rate
% - Param.nf = FFT n. of points
% - Param.av = averaging type: 1=arithm.mean, 2=median
% - Param.fn = n. of FFTs in each snapshot
% - Param.of = offset in samples between successive FFTs
% - Param.zp = zero pad samples added for last FFT
% - Param.sn = current snapshot number
% - Param.re = pressure reference (normally 1 microPa)
% - Param.vb = verbose mode (0=silent, 1=verbose)
% - Param.pp = pause after each plot in s
% - Param.figwav = figure handle for WAV plot, used in debug mode
% - Param.figpsdb = figure handle for broadband PSD plot
% - Param.figpsdd = figure handle for decidecade PSD plot
% - Param.procpsdb = broadband PSD processing selected (boolean)
% - Param.procspld = decidecade SPL processing selected (boolean)

fftno = ceil(Param.fn); % n. of FFTs, including fractional last one
PREF = Param.re;

VERB = Param.vb;
PLTPAUSE = Param.pp;
PROCPSDB = Param.procpsdb; % broadband PSD
PROCSPLD = Param.procspld; % decidecade SPL
PROCSPLB = Param.procsplb; % broadband SPL

% plot titles and labels
CURTIT = 'Snap %d, FFT %d - high/current/low'; % 2 arguments
AVGTIT = 'Snap %d, total %d FFTs - high/average/low'; % 2 arguments
BBXLAB = 'Frequency / Hz';
BBYLAB = 'PSD / dB re 1 \muPa^2/Hz';
DDXLAB = 'Center frequency / Hz';
DDYLAB = 'Decidecade SPL / dB re 1 \muPa^2';

figWAV = Param.figwav; % WAV figure handle
if PROCPSDB
    figPSDB = Param.figpsdb; % broadband PSD figure handle
end
if PROCSPLD
    figSPLD = Param.figspld; % decidecade PSD figure handle
end

DOTN = 10; % display '.' character every N cycles to show progress

fbinnorm = (Param.sr/Param.nf); % normalize for f.bin width (linear)
zeropadnorm = (1 - Param.zp/Param.nf); % normalize for zero pad: UNUSED

if VERB > 0    % add details if verbose mode on
    fprintf('\nParam: sr=%g, nf=%d, fn=%g, of=%g, zp=%d', ...
        Param.sr, ...
        Param.nf, ...
        Param.fn, ...
        Param.of, ...
        Param.zp)
    fprintf('\nNorm: fbin=%g, zp=%g',fbinnorm,zeropadnorm); 
end

bbpxxmat = zeros(Param.nf/2 + 1, fftno); % allocate BB buffer matrix

if PROCSPLD
    ddf = DDmat(:,1); % vector of DD center freq.s, used in plots
    ddpxxmat = zeros(size(DDmat,1), fftno); % allocate DD buffer
end

%%%%%%%%%%%%%%%%% loop over n. of FFTs within snapshot %%%%%%%%%%%%%%%%
for i = 1:fftno % last portion might be fractional
    
    startSmp = 1 + (i-1) * Param.of;
    stopSmp  = startSmp + Param.nf - 1;
    
    if stopSmp > length(Snap.Y) % last FFT vector may exceed snapshot:
        % slide back to take full last portion of snapshot
        lastsubint = Snap.Y(startSmp:length(Snap.Y)); % fewer samples
        zeropadint = zeros(Param.zp, 1); % zeropad, row-wise
%        Y = Snap.Y(length(Snap.Y) - Param.nf + 1:length(Snap.Y));
        Ynoslide = vertcat(lastsubint, zeropadint); % zero padded subint:
                                            % only for visual purpose
        ovlap = stopSmp - length(Snap.Y);
        startSmp = length(Snap.Y) - Param.nf + 1; % slide back % update
        stopSmp  = length(Snap.Y); % update to end of snapshot
        zpString = sprintf(" len=%d, ovl=%d, ",length(lastsubint),ovlap);
        islastsubint = true;
    else
        zpString = '';
        islastsubint = false;
    end % %%%%%%%%%%%%%%%%%% last subinterval

    Y = Snap.Y(startSmp:stopSmp); %%%%%%%%%%%% build FFT subinterval

    %%%%%%%%%%%%%% DEBUG: draw Y (zeropad and slide back)
    if VERB == 2
        % figure(figWAV);
        set(0,'CurrentFigure',figWAV); % draw on waveform figure
        if islastsubint
            plot(1:Param.nf,Ynoslide) % show zero padded subint
            fprintf(' orig+zeropad...')
            pause % hit any key
            fprintf(' slide back...')
        end
        plot(1:Param.nf,Y)
        pause % hit any key
    end
    %%%%%%%%%%%%%% END DEBUG
    
    if VERB > 0 % display progress info
        fprintf('\n%d: [%d %d] %s',i,startSmp,stopSmp,zpString); 
    else
        if mod(i,DOTN) == 0
            fprintf('.'); % silent progress, only print dot every DOTN
        end
    end

    %%%%%%%%%%%%% compute raw PSD through periodogram on windowed data;
    % output is normalized to power of window specified as parameter

    [pxx,bbf] = periodogram(Y, Win, Param.nf, Param.sr);
    
    % no need to normalize: subinterval length is fixed using slideback
    % if i == fftno && Param.zp > 0 % last FFT with zero padding
    %     pxx = pxx / zeropadnorm; % normalize to zeropad length
    % end

    bbpxxSF = (pxx .* (BBSF'.^2)); % apply BB calibration scale factor

    bbpxxmat(:,i) = bbpxxSF; % update column vector in buffer matrix

    %%%%%%%% BB min-max bounds, to be repeated over DDecs if enabled
    bbpxxdb = 10*log10(bbpxxSF) - 20*log10(PREF); % convert to dB
    if VERB > 0
        fprintf(' max = %g dB',max(bbpxxdb)); % disp. max if verbose
    end

    if i == 1   % initialize max & min
        bbpxxdbhi = bbpxxdb;
        bbpxxdblo = bbpxxdb;
    else           % update max & min
        for n = 1 : length(bbpxxdb)  
            if bbpxxdb(n) > bbpxxdbhi(n)
                bbpxxdbhi(n) = bbpxxdb(n);
            end
            if bbpxxdb(n) < bbpxxdblo(n)
                bbpxxdblo(n) = bbpxxdb(n);
            end
        end
    end % BB min-max bounds

    if PROCSPLD %%%%%%%%% decidecade min-max bounds, repeat from above.
        % For each DD band, add BB elements & divide by DD bandwidth.
        for ddidx = 1:size(DDmat,1) % DD index
            fvil = int32(DDmat(ddidx,6)); % lower DD limit f.vect.idx
            fvih = int32(DDmat(ddidx,7)); % upper DD limit f.vect.idx
            ddlow = DDmat(ddidx,2);       % nominal DD lower limit in Hz
            ddhi  = DDmat(ddidx,3);       % nominal DD upper limit in Hz
            
            % f bin normalization factor = s.rate/n.fft: if = 1
            % then each element is exactly 1 Hz wide, otherwise need to
            % account for bin width to compute total power.
            
            %  OLD v. 3.x: decidecade power was averaged over DD band;
            %  removed division by bandwidth as of v. 4.x
%           ddpxx(ddidx,i) = sum(pxx(fvil:fvih))*fbinnorm/(ddhi-ddlow);
            ddpxx(ddidx,i) = sum(pxx(fvil:fvih)) * fbinnorm;

            % OLD     ddpxx(ddidx,i) = mean(pxx(fvil:fvih));
            % in the old version above, taking the mean was equivalent to
            % summing and dividing by (b.width in Hz / bin width in Hz)
            % where the quantity in () is the n. of DD band elements
        end

        ddpxxSF(:,i) = ddpxx(:,i) .* (DDSF.^2); % decidecade calib.SF
        ddpxxdb = 10*log10(ddpxxSF(:,i)) - 20*log10(PREF); % to dB

        if i == 1   % initialize max & min on first FFT cycle
            ddpxxdbhi = ddpxxdb;
            ddpxxdblo = ddpxxdb;
        else
            for n = 1 : size(ddpxxmat,1)   % update max & min
                if ddpxxdb(n) > ddpxxdbhi(n)
                    ddpxxdbhi(n) = ddpxxdb(n);
                end
                if ddpxxdb(n) < ddpxxdblo(n)
                    ddpxxdblo(n) = ddpxxdb(n);
                end
            end
        end
    end % decidecade min-max bounds

    %%%%%% plot hi-current-low only on enabled figures,
         % although BB data are always present
    if PROCPSDB
        % figure(figPSDB);
        set(0,'CurrentFigure',figPSDB); % draw on BB PSD figure
        semilogx(bbf,bbpxxdbhi,'r',bbf,bbpxxdblo,'g',bbf,bbpxxdb,'b' )
        grid on
        title(sprintf(CURTIT,Param.sn,i));
        xlabel(BBXLAB)
        ylabel(BBYLAB)
    end
    if PROCSPLD
        % figure(figSPLD);
        set(0,'CurrentFigure',figSPLD); % draw on DD PSD figure
        semilogx(ddf,ddpxxdbhi,':or', ...
                 ddf,ddpxxdblo,':og', ...
                 ddf,ddpxxdb,'ob' )
        grid on
        title(sprintf(CURTIT,Param.sn,i));
        xlabel(DDXLAB)
        ylabel(DDYLAB)
    end
    
    if VERB == 2
        pause %%%%%% DEBUG: hit any key to continue
    else
        pause(PLTPAUSE/10) % normal run should be fast, short pause
    end

end %%%%%%%%%%%%%% end FFT loop %%%%%%%%%%%%%%%%%

if PROCSPLB || PROCPSDB   % broadband SPL or PSD
    switch Param.av % return average on FFT subint. linear data
        case 1 % arithmetic mean
            bbpxxavg = mean(bbpxxmat,2); % V.3.x BUG: was on bbpxxSF
        case 2 % median
            bbpxxavg = median(bbpxxmat,2); % as above
    end
end

if PROCSPLB     % From v. 4.0:  return full spectrum power as
    % sum of components weighted by f-dependent SF, instead of 
    % using wave rms and mean SF as in earlier main v. 5.xx
    fbpxxavg = sum(bbpxxavg,1) * fbinnorm;
else
    fbpxxavg = []; % return void if not selected in main, see below.
end

%%%%%% plot hi-low-averages for BB, DD, if enabled

if PROCPSDB % broadband

    bbpxxdbavg = 10*log10(bbpxxavg) - 20*log10(PREF); 
    figure(figPSDB); % put figure on top before return
    % set(0,'CurrentFigure',figPSDB);
    semilogx(bbf,bbpxxdbhi,'r', bbf,bbpxxdblo,'g',bbf,bbpxxdbavg,'b')
    grid on
    title(sprintf(AVGTIT,Param.sn,i));
    xlabel(BBXLAB)
    ylabel(BBYLAB)
else
    bbpxxavg = []; % return void if BB option not selected in main:
    % need to define a return value even if option is not selected,
    % otherwise an error is generated in main.
end

if PROCSPLD
    switch Param.av % return average on linear data
        case 1 % arithmetic mean
            ddpxxavg =   mean(ddpxxSF,2);
        case 2 % median
            ddpxxavg = median(ddpxxSF,2);
    end
    ddpxxdbavg = 10*log10(ddpxxavg) - 20*log10(PREF); 
    figure(figSPLD); % put figure on top before return
    % set(0,'CurrentFigure',figSPLD);
    semilogx(ddf,ddpxxdbhi,':or',ddf,ddpxxdblo,':og') % plot hi-low
    hold on
    semilogx(ddf,ddpxxdbavg,'ob','MarkerFaceColor','b') % plot avg
    hold off
    grid on
    title(sprintf(AVGTIT,Param.sn,i));
    xlabel(DDXLAB)
    ylabel(DDYLAB)
else
    ddpxxavg = []; % return void if DD option not selected (see BB above)
end

if VERB == 2
    pause % DEBUG: hit any key to continue
else
    pause(PLTPAUSE) % normal run, pause as user specified
end
return

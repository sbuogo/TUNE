% T U N E = Tool for Underwater Noise Evaluation: main function called by 
%           graphical user interface
% 
% DESCRIPTION
% Function is called with parameter setup as argument, and returns
% exit status, error message, and figure handles to user interface.
% Setup also specifies a directory with all audio WAV files to be
% processed, another one to store plots and results, and a calibration
% file in text format with row-wise pairs Hz-dB or Hz-SF (scale factor).
% Processing options: sound pressure level (SPL: broadband, decidecade,
% zero-to-peak); power spectral density (PSD). Processing is done in 
% intervals (snapshots) of given duration between given start and end 
% times, or up to end of file. Snapshots may be consecutive, overlapping 
% or separated by a given time. DC component is removed and high/low/
% bandpass filter may be applied before processing. For each snapshot,
% FFT is computed in N-sample overlapping subintervals (with N power of 2) 
% and averaged over all subintervals to give a snapshot spectrum. 
% All snapshot spectra are then averaged (arithmetic mean or median) and 
% plotted in dB units with optional high-low percentiles.
% Decidecade analysis is done in a set of consecutive bands by specifying
% start and end nominal center frequencies (ISO preferred frequencies).
% Decidecade limits are computed using formula 10^(m/20) with
% m = [2n-1, 2n+1] for each n-index center frequency (index 0 = 1 Hz).
% Results are in dB relative to 1 microPa standard reference pressure.
% Info may be logged to file during processing. Output plots are shown
% during processing for each selected option, and final ones may be
% saved to graphic files. Data may also be saved to text files.
% Setup parameters may be saved to, or loaded from, a given *.mat file.
%
% The following files must be in the same directory as this script:
%  - TUNE[xxx].mlapp	 where [xxx] is current version number
%  - TUNEsub[nn].m       where [nn] is required version (see below)
%
% NOTE ON MATLAB VERSIONS
% Present version requires Matlab R2012 or later. Formatted text output
% is done using dlmwrite() function up to R2018, otherwise the 
% recommended function writematrix(), introduced in R2019, is used.
% 
% REVISION HISTORY
% Up to v.3.2: development versions
% v.3.3: Correct handling of first&last TOB, no illegal array index.
% v.4.1: Batch mode, param.s & opt.s declared here as static.
% v.4.2: Imported Knudsen curve overlay from 'MSPSD21.m'.
% v.4.3: Use correct file path separator for Win or Linux ('\' or '/'),
%        added option for simple/box TOB plot and max,avg in plot legend.
%        Renamed Y axis to '3rd-octave SPL' in TOB plot (was same as PSD).
%        New parameters passed to subroutines: p_ref and verb.mode.
% v.4.4: Results averaging now made on linear values, not on dB values
%        (requires broadbandspl() v.2.3, tobspl() v.2.4 or later).
% v.4.5: Changed function 'contains' to 'strcmp' due to incompatibility
%        with Matlab versions prior to R2016.
% v.4.6: Improved broadband plot with mean/median and percentiles.
%        Fixed start folder in file selection windows, and minor bugs.
% v.4.7: Revised terminology, renamed function 'broadbandspl' to 'evalpsd'
%        and 'tobspl' to 'evaltob'. Changed related filenames accordingly.
%        Added warning to increase FFT points if unable to resolve TOB.
% v.4.8: option for SPL, use mean cal.SF over full band, wave plot in Pa
%        optional highpass filter; 60 dB stopb.atten, 0.1 dB ripple, 
%        trans.band = 0.15*fc. Multiple processing types now possible, 
%        options are single boolean. Waves are always displayed first
%        in Pa using calibration, each plot is on separate figures.
% v5.02: Changed to function called by GUI developed with App Designer.
%        Parameters are passed using setup (struct type) as argument.
%        Added optional highpass filter before processing.
%        P_ref (1 microPa) no longer user defined, now static value.
%        Added command output log (diary). Return value: 0=ok, 1=error,
%        return exit status string on error. Plot titles now optional.
%        Setup optionally logged to file.
% v5.03: Plot labels, titles and legends defined globally.
% v5.04: Figures positioned near screen bottom left corner, height
%        reduced if exceeding over screen top.
% v5.05: Added dialog box for info display. Added set(currentfigure)
%        upon each plot to avoid focus stealing when clicking on another.
%        Subroutines 'evalPSd','evalTOB' changed, added version number.
%        Requires subroutines 'evalPSD' v.2.6, 'evalTOB' v.2.7.
%        Handles to current figures passed to GUI to delete them on exit.
% v5.06: Merged subroutines 'evalPSD' and 'evalTOB' to 'evalPSDBBDD' v.3.0
%        replacing one-third octave with decidecade f intervals.
%        Accepts multichannel audio files processing one channel at a time.
% v5.10: Setup may be saved to *.mat file or loaded from existing one.
%        Added optional lowpass filter along with existing highpass
%        using same parameters: use bandpass mode if both are selected.
% v5.11: Fixed error in subroutine 'eval*' due to bad return value.
%        Last FFT subinterval no longer handled using zeropad, instead 
%        slide back to last full snap subinterval. Use subroutine v.3.2.
% v5.12: Remember past input/output folders and cal. file on new session.
% v5.13: Name change to TUNE, removed zeropad normalization in subroutine
%        v.3.3. Removed backwards compatibility prior to Matlab R2012.
%        Added header lines to output text files. Renamed DD-PSD to DD-SPL
%        and related parameters accordingly. New subroutine v.3.3. Old
%        setup files are no longer supported (no check is done).
% v5.14: Fixed bug with decimal comma in calibration file, added warning
%        and figure plot. Exclusive snap skip/overlap passed from GUI.
%        Version n. declared in GUI to avoid using other setup versions.
% v6.01: Redefined broadband and decidecade SPL calculations: broadband
%	 now computed summing all PSD components wheighted by individual SF
%	 (was done by taking wave rms and mean SF). Decidecade SPL is now
%	 the sum of components within each band, and not its average (was
%	 done by summing and dividing by bandwidth, now it is only a sum).
% v6.02: Fixed bug on cal.file interpolation, now linear between f points
%    as 'nearest' method could not be used for extrapolation only.
% v6.11: Added option for zero-to-peak level. Setup file from a different
%    version can now be selected, a GUI warning is issued if not valid.
%    Subroutine and m-file renamed to "TUNEsubNN" with NN version number.
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
% VERSION 6.11 | 03-oct-2025 

%%%%%%% NOTE: upon version change, update function name.

function [retVal,exStat,figHandles] = TUNE611main(Setup)

% VERS = 'X.XX'  since v. 5.14 it is included in setup 
PROGNAME = 'TUNE';

% setup argument contains all parameters defined by user in GUI
% NOTE: new definitions as of v. 5.14, no backwards compatibility
VERS        = Setup.PROGVERS;
PROCSPLB    = Setup.PROCSPLB; % broadband SPL (was PROCSPL)
PROCPSDB    = Setup.PROCPSDB; % broadband PSD
PROCSPLD    = Setup.PROCSPLD; % decidecade SPL (was PROCPSDD)
PROCSPLP    = Setup.PROCSPLP; % zero-to-peak SPL (new in v. 6.11)
BOXPLT      = Setup.BOXPLT;
STARTTIME   = Setup.STARTTIME;
STOPTIME    = Setup.STOPTIME;
SNAPLEN     = Setup.SNAPLEN;
SNAPOVLAP   = Setup.SNAPOVLAP;
SNAPSKIP    = Setup.SNAPSKIP;
NFFT        = Setup.NFFT;
FFTOVLAP    = Setup.FFTOVLAP;
WTYPE       = Setup.WTYPE;
AVGTYPE     = Setup.AVGTYPE;
HIPASSF     = Setup.HIPASSF;
LOPASSF     = Setup.LOPASSF;
PRCTILEL    = Setup.PRCTILEL;
PRCTILEH    = Setup.PRCTILEH;
DDIDXL      = Setup.DDIDXL;
DDIDXH      = Setup.DDIDXH;
DRAWALL     = Setup.DRAWALL;
SAVEJPG     = Setup.SAVEJPG;
SAVETXT     = Setup.SAVETXT;
XAXAUTO     = Setup.XAXAUTO;
YAXAUTO     = Setup.YAXAUTO;
XAXMIN      = Setup.XAXMIN;
XAXMAX      = Setup.XAXMAX;
YAXMIN      = Setup.YAXMIN;
YAXMAX      = Setup.YAXMAX;
LEGND       = Setup.LEGND;
FIGSIZE     = Setup.FIGSIZE;
VERB        = Setup.VERB;
PLTPAUSE    = Setup.PLTPAUSE;
KNDSDRAW    = Setup.KNDSDRAW;
KNDSSIDX    = Setup.KNDSSIDX;
% next added in v. 5.02
INPUTPATH   = Setup.INPATH; 
OUTPUTPATH  = Setup.OUTPATH;
CALFILEPATH = Setup.CALFPATH;
LOGTOFILE   = Setup.LOGFILE; % save command output to file
TITL        = Setup.TITL; % show title in all figures
LOGSETUP    = Setup.LOGSETUP; % add setup in log file
% next added in v. 5.05
FREEFIG     = Setup.NORESIZ; % keep current figure size and positions
% end of setup argument fields

%%% Decidecade index definitions, kept here as a reference.
%%% NOTE: center frequencies are nominal ISO, while band limits are 
%%% computed using power-of-ten formula.
%
%   1=1Hz     11=10Hz    21=100Hz   31=1000Hz   41=10kHz    51=100kHz
%   2=1.25 Hz 12=12.5Hz  22=125Hz   32=1250Hz   42=12.5kHz  52=125kHz
%   ...       13=16Hz    23=160Hz   33=1600Hz   43=16kHz    ...
%   ...       14=20Hz    24=200Hz   34=2000Hz   44=20kHz    ...
%   ...       15=25Hz    25=250Hz   35=2500Hz   45=25kHz
%             16=31.5Hz  26=315Hz   36=3150Hz   46=31.5kHz
%             17=40Hz    27=400Hz   37=4000Hz   47=40kHz
%             18=50Hz    28=500Hz   38=5000Hz   48=50kHz
%   ...       19=63Hz    29=630Hz   39=6300Hz   49=63kHz    ...
%   10=8 Hz   20=80Hz    30=800Hz   40=8000Hz   50=80kHz    60=800kHz

%%%%%% Figure labels, titles and legend - may be changed if requested
FigLab.XLABWAV  = 'Time / s';
FigLab.XLABSPLB = FigLab.XLABWAV;
FigLab.XLABSPLP = FigLab.XLABWAV;
FigLab.XLABPSDB = 'Frequency / Hz';
FigLab.XLABSPLD = 'Center frequency / Hz';
FigLab.YLABWAV  = 'Sound pressure / Pa';
FigLab.YLABSPLB = 'Sound pressure level, L_{p,rms} / dB re 1 \muPa^2';
FigLab.YLABSPLP = 'Peak sound pressure level, L_{p,pk} / dB re 1 \muPa';
FigLab.YLABPSDB = 'Power spectral density, L_{p,f} / dB re 1 \muPa^2/Hz';
FigLab.YLABSPLD = 'Decidecade SPL, L_{p,dd} / dB re 1 \muPa^2';
% titles are passed to sprintf() with placeholders for arguments:
% observe n. of placeholders for arguments: wav=3;spl-bb=4;psd,spl-dd=5
FigLab.TITLWAV  = '%s: snap %d/%d';
FigLab.TITLSPLB = '%s: %.0f-%.0f s, %g s snap';
FigLab.TITLSPLP = FigLab.TITLSPLB;
FigLab.TITLPSDB = '%s: %.0f-%.0f s, %g s snap, %d pts FFT';
FigLab.TITLSPLD = FigLab.TITLPSDB;
FigLab.LGNDPSDB = [strcat(num2str(PRCTILEL),"th percentile"); ...
                   strcat(num2str(PRCTILEH),"th percentile") ];
FigLab.LGNDSPLD = FigLab.LGNDPSDB;
% note: all strings in boxplot legend must be equal length
FigLab.LGNDBOX = ['Box edges = 25th-75th percentiles'; ...
                  'Box center line = median         '; ...
                  'Dashed vertical line = min-max   '; ...
                  'Cross = outlier                  ' ];

%%%%%% static global parameters - do not change
PREF      = 1e-6;   % Reference pressure in Pa (standard: 1 microPa)
ISOD = [1 1.25 1.6 2 2.5 3.15 4 5 6.3 8]; % ISO preferred f in a decade
ISOF = [ISOD ISOD*10 ISOD*100 ISOD*1e3 ISOD*1e4 ISOD*1e5]; % 1Hz - 800kHz
% Knudsen's reference noise levels for sea states 0 to 6 (Ref: Lurton)
KNDLEVS = [44.5  55  61.5  64.5  66.5  68.5  70]; % @1kHz, SS0 to SS6
SSLABELS = ['SS0';'SS1';'SS2';'SS3';'SS4';'SS5';'SS6']; % Knudsen labels
KNDDEC = (-17);   % NL = KNDLEVS + KNDDEC*log10(f in kHz) for f > KNDCUT
KNDCUT = 1e3;     % level assumed constant up to this frequency in Hz
TXTXPOS = 1.5*XAXMIN;   % X position for knudsen curve labels
TXTYOFFS = 1.5;         % label Y offset in dB with respect to line
FIGWAV = 1; % figure handle numbers 
FIGSPLB = 2; 
FIGPSDB = 3; 
FIGSPLD = 4; 
FIGSPLP = 5; 
% check for any reserved n. in GUI before adding more
%FIGINFO = 9; % additional figure for info display (changed to msgbox)
DISPINFO = true; % true=show info msgbox figure, false=do not show
FIGXYPOS = [860 130]; % initial X Y fig.pos., pix from bottom left corner
FIGXYOFFS = [50 -35]; % incremental X Y offset for each figure in pix
FIGINFOXY = [75 545]; % 'info' modal figure position, see GUI properties
FILTBOXPOS = [0.3 0.3 0.4 0.1]; % "filtering..." box msg. pos.: X Y W H
FIGDRV = '-djpeg'; % graphics driver for figure output (e.g. jpeg)
FIGRES = '-r300'; % output figure resolution in dpi
TXTDECS = 3; % n. of decimals in output text data (default = 13)
PEAKFMIN = 10; % min f in Hz to detect spectrum peak (for peak SPL)
LOGFILENAME = [PROGNAME '_log.txt']; % log file to save command output,
                              % must match with that in GUI code
HDR    = ['%% ' PROGNAME '  ' VERS '\n']; % output text file header
HDRSPLB = [HDR '%% Broadband SPL \n%% t/s \tSPL/dB\n'];
HDRSPLP = [HDR '%% Zero-to-peak SPL \n%% t/s \tSPL/dB\n'];
HDRPSDB = [HDR '%% Broadband PSD \n%% f/Hz\t[Avg\tPct-\tPct+] / dB\n'];
HDRSPLD = [HDR '%% Decidecade SPL\n%% f/Hz\t[Avg\tPct-\tPct+] / dB\n'];
%%%%%% end of static global parameters 

if ispc  % file path separator: '\' for Windows, '/' for Linux
    PathSeparator = '\';
else
    PathSeparator = '/';
end
matlabVers = version('-release');
versYear = str2num( matlabVers(1:4) ); % skip trailing 'a' or 'b'

retVal = 0; % clear return value to 0 (no error) at beginning
exStat = ''; % on error, return error description string
figHandles = ''; % do not return figure handles on error

if versYear < 2012
    exStat = 'Matlab version R2012 or newer required. Exiting.';
    disp(exStat)
    retVal = 1;
    msgbox(exStat,'Error','error','modal');
    return
end

timeOfStart = datetime('now');          
if LOGTOFILE
    diary([OUTPUTPATH,PathSeparator,LOGFILENAME])
    disp(char(61*ones(1,80))) % separator line: ascii 61 = '='
    fprintf('Start: %s\n',timeOfStart)
else
    diary off
end

% display header in command window and logfile if enabled
hdr=['| T U N E  =  Tool for Underwater Noise Evaluation  ',VERS,' |'];
hdrb = char([111 45*ones(1,length(hdr)-2) 111]); % ascii 45='-', 111='o' 
disp([hdrb; hdr; hdrb]);

%%%%% select folders input, output, and calib.file
InputFolder = INPUTPATH;
OutputFolder = OUTPUTPATH;
wavPath = [InputFolder PathSeparator];  
outPath = [OutputFolder PathSeparator];
wFiles = dir(InputFolder);  

validWFiles = 0;
for idx = 3:length(wFiles)   % count all valid WAV files skipping current
                             % and parent directory which come first
    len = length(wFiles(idx).name);
    if len > 3
        ext = lower(wFiles(idx).name(len-3:len)); % .wav for valid file
    else
        ext = [];   % not a valid filename with 3-char extension
    end
    if wFiles(idx).isdir == 0 && strcmp(ext,'.wav')
        validWFiles = validWFiles+1;
    end
end
fprintf('Input folder %s:\n  found %d WAV files out of %d total\n',...
            InputFolder,validWFiles,length(wFiles)-2);
fprintf('Output folder %s\n', OutputFolder);
        
%%%%% load calibration file
CalFilePath = CALFILEPATH;
if CalFilePath == 0 % user canceled
    exStat = 'Error: must select a valid calibration file.';
    disp(exStat)
    retVal = 1;
    msgbox(exStat,'Error','error','modal');
    return
else
    if ~isfile(CalFilePath) % initial condition, no file
        exStat = 'Error: missing calibration file.';
        disp(exStat)
        retVal = 1;
        msgbox(exStat,'Error','error','modal');
        return
    end
end
[~,CalFileName,CalFileExt] = fileparts(CalFilePath); % folder is unused
calTab = readtable(CalFilePath); % blank lines give NaN
calMat = calTab{:,:}; % convert table into matrix for further processing
calMat = rmmissing(calMat); % remove NaN if blank lines in table
calPts = size(calMat,1);
if calPts < 3
    exStat = 'Error: at least 3 points required in calibration file.';
    disp(exStat)
    retVal = 1;
    msgbox(exStat,'Error','error','modal');
    return
end

fprintf('Calibration file %s:\n  %d points spanning [%g - %g] Hz, ',...
    CalFilePath, calPts, calMat(1,1), calMat(calPts,1) );

if max(calMat(:,2)) > 0 % file is in SF(Pa) units, convert to dB
    fprintf('SF(Pa) units\n');
    calMat(:,2) = 20*log10( 1./calMat(:,2) ) - 120;
else
    fprintf('dB units\n'); % leave as is
end

refdb = 20 * log10(PREF); % all levels are relative to reference p

%%%%%%%%%%%%%%%%%%%%% main loop over audio files %%%%%%%%%%%%%%%%%%%%%%%%

for wFileN = 3:length(wFiles)

    len = length(wFiles(wFileN).name);
    if len > 3
        ext = lower(wFiles(wFileN).name(len-3:len)); % .wav for valid file
    else
        ext = [];   % no filename extension
    end

    if wFiles(wFileN).isdir > 0 || ~strcmp(ext,'.wav')
            fprintf('\n"%s": not a WAV file\n',wFiles(wFileN).name)
        continue    % skip to next audio file
    else
        % valid file, proceed with processing
        WavFileName = wFiles(wFileN).name;
    end
    fprintf('\nInput file %d of %d: "%s"\n',...
                        wFileN-2, length(wFiles)-2, WavFileName);

    % initialize figures, later on to be called by set(0,'currentfigure')
    figWav=figure(FIGWAV); set(gcf,'Name','Waveform');
            set(gcf,'NumberTitle','off');
    if PROCSPLB
        figSplB=figure(FIGSPLB); set(gcf,'Name','Broadband SPL');
            set(gcf,'NumberTitle','off');
    end
    if PROCPSDB
        figPsdB=figure(FIGPSDB); set(gcf,'Name','Broadband PSD');
            set(gcf,'NumberTitle','off');
    end
    if PROCSPLD
        figSplD=figure(FIGSPLD); set(gcf,'Name','Decidecade SPL');
            set(gcf,'NumberTitle','off');
    end
    if PROCSPLP
        figSplP=figure(FIGSPLP); set(gcf,'Name','Zero-to-peak SPL');
            set(gcf,'NumberTitle','off');
    end

    WavFile = [wavPath WavFileName]; % not affected by trailing '(chN)'
                                     % if multichannel

    [F.path,F.name,F.ext] = fileparts(WavFile);
    FileNameNoExt = F.name; % output files share same filename

    WavFileInfo = audioinfo(WavFile); % contains file header info

    smpRate = WavFileInfo.SampleRate;
    nSmp    = WavFileInfo.TotalSamples;
    fileDur = WavFileInfo.Duration; % in seconds
    fileDurStr = duration(0,0,fileDur,'Format','hh:mm:ss.SSS');
    nChans  = WavFileInfo.NumChannels;
    chanStr = 'channel';
    if nChans > 1
        chanStr = [chanStr 's'];
    end
    nBits = WavFileInfo.BitsPerSample;
    if nSmp > 0
        fprintf(' Rate: %g kHz, %.3f Msamples,', smpRate/1e3, nSmp/1e6);
        fprintf(' duration: %s (%.3f s)\n', fileDurStr, fileDur);
        fprintf(' %d bits per sample, %d %s\n', nBits, nChans, chanStr);
        if not(isinf(STOPTIME)) && STOPTIME > fileDur
            exStat = '* Warning: stop time exceeds file duration.';
            disp(exStat)
            retVal = 1;
            validWFiles = validWFiles - 1; % files actually processed
            continue % does not exit program, skip to next file
        end
    else
        exStat = '* Warning: acoustic file seems to be void.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % does not exit program, skip to next file
    end

    WavFileNameOrig = WavFileName; % to append (chN) if multichannel
    FileNameNoExtOrig = FileNameNoExt;

    if HIPASSF >= smpRate || LOPASSF >= smpRate
        exStat = '* Warning: filter cut f exceeds sample rate.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % does not exit program, skip to next file
    end

    
    %%%%%%%%%%%%%%%%% loop over channels if more than 1 %%%%%%%%%%%%%%%

    for curChan = 1:nChans 
        % loop index also used to avoid repeating global info display

    if nChans > 1 % stereo or multichannel, append (chN) to filename
        % text message output done later on, after global setup
        WavFileName   = [ WavFileNameOrig   sprintf('(ch%d)',curChan)];
        FileNameNoExt = [ FileNameNoExtOrig sprintf('(ch%d)',curChan)];
    end

    if isinf(STOPTIME)
        stopTime = fileDur;
    else
        stopTime  = STOPTIME;
    end
    snapLen   = SNAPLEN;
    snapOvlap = SNAPOVLAP;
    snapSkip  = SNAPSKIP;

    if snapLen > fileDur
        exStat = '* Warning: snapshot exceeds file duration.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % does not exit program, skip to next file
    end

    % manage snapshot overlap and skip time
    nSnapOvlap = (100/(100-snapOvlap)); % snap.s starting within snapLen
    % nSnap is total no. of snapshots accounting for possible overlap
    nSnapExact = nSnapOvlap*((stopTime-STARTTIME)/snapLen - 1);
    nSnap = floor(nSnapExact) + 1;
    snapStep = snapLen / nSnapOvlap; % time step between snapshots
    % override if skip > 0 (get rid of overlap, should be 0 anyway)
    if snapSkip > 0
        nSnap = floor((stopTime-STARTTIME)/(snapLen+snapSkip));
        snapStep = snapLen + snapSkip;
    end

    % build matrix of [start,stop] samples for each snap
    snapSmp  = zeros(nSnap,2); % 1st col=start, 2nd col=stop
    for idx = 1:nSnap % set start & stop smp for each snap
        snapSmp(idx,1) = round((STARTTIME+(idx-1)*snapStep)*smpRate+1);
        snapSmp(idx,2) = round(snapSmp(idx,1) + snapLen*smpRate - 1);
    end
    snapTime = snapSmp/smpRate; % start & stop time for each snap

    % define snapshot timebase
    Snap.T = 0 : 1/smpRate : snapLen-(1/smpRate); % struct: Snap

    if curChan == 1 % do not repeat info display for next channels
        fprintf(' Timebase: [%.1f-%.1f] s, %d snapshots, %.1f s each', ...
            STARTTIME,stopTime,nSnap,snapLen);
        if SNAPOVLAP > 0
            fprintf(', %d%% overlap', SNAPOVLAP);
        end
        if SNAPSKIP > 0
            fprintf(', separated by %.1f s', SNAPSKIP);
        end
        fprintf('\n');
    end

    % manage FFT overlap
    nfftOvlap = (100/(100-FFTOVLAP)); % FFTs starting within each nFFT
    % n. of FFTs in snapshot, accounting for possible overlap:
    snapFFTnexact = nfftOvlap*(snapLen*smpRate/NFFT - 1) + 1;
    snapFFTn = ceil(snapFFTnexact); % round up for zeropad
    offsFFT = NFFT / nfftOvlap; % incremental offset of each FFT
    % total samples spanned by FFT, accounting for possible overlap:
    spanFFT = NFFT + (snapFFTn-1)*offsFFT; 

    if snapFFTnexact < 1    % cannot manage FFT if points are more than
                            % samples in snap: slideback would fail
        exStat = ...
        '* Warning: FFT does not fit into snapshot: reduce FFT points.';
        disp(exStat)
        retVal = 1;
        validWFiles = validWFiles - 1; % files actually processed
        continue % do not exit program, try next file
    end

    if snapFFTnexact-fix(snapFFTnexact) > 0 % noninteger FFTn: zero pad
        zeroPad = spanFFT - snapLen*smpRate;
    else   % exact n. of FFTs in snapshot: no zero pad
        zeroPad = 0;  
    end

    % broadband f vector ranges from f=0 to fmax=(smp.rate)/2
    fVect = linspace(0, smpRate/2, NFFT/2 + 1);

    % Interpolate calibration data (in dB) in FFT frequency points.
    % No longer use 'nearest' method with 'extrap', use 'linear'
    % then extrapolate forcing fixed values on both sides.
    % NOTE: to be done for each file as it depends on sample rate
    
    CaldB = interp1(calMat(:,1),calMat(:,2),fVect,'linear');
    fLowEnd  = fVect < calMat(1,1); % logical arrays for extrapolation
    fHighEnd = fVect > calMat(size(calMat,1),1);
    CaldB(fLowEnd)  = calMat(1,2); % extrapolate first value...
    CaldB(fHighEnd) = calMat(size(calMat,1),2); % and last value
    CalSF  = 1./(10.^((CaldB+120)/20)); % convert to Scale Factor in Pa
    meanCalSF = mean(CalSF);  % for waveform plot

    wType = WTYPE;
    switch wType
        case 0
            FFTW = ones(NFFT,1); % rectangular (no window)
        case 1
            FFTW = hamming(NFFT);
        case 2
            FFTW = blackman(NFFT);
        case 3
            FFTW = hann(NFFT);
        otherwise
            FFTW = ones(NFFT,1);
            if VERB > 0
                disp('* Invalid window, using rectangular');
            end
    end
    avgType = AVGTYPE;  % arithmetic mean or median

    % set of parameters to be passed on to external functions
    Param.sr = smpRate;
    Param.nf = NFFT;    % FFT n. of points 
    Param.av = avgType;
    Param.fn = snapFFTnexact; % n. of FFTs in each snap
    Param.of = offsFFT;  % offset between FFTs
    Param.zp = zeroPad;  % zero padding for last FFT
    Param.re = PREF;    % reference pressure
    Param.vb = VERB;    % verbose mode
    Param.pp = PLTPAUSE; % pause after each plot
    Param.procpsdb = PROCPSDB; % boolean
    Param.procspld = PROCSPLD; % boolean
    Param.procsplb = PROCSPLB; % boolean, added in v. 6.01
    Param.procsplp = PROCSPLP; % boolean, added in v. 6.11
    Param.peakfmin = PEAKFMIN; % in Hz, for peakSPL option (v. 6.11)
    Param.figwav = figWav; % wav figure handle
    if PROCPSDB
        Param.figpsdb = figPsdB; % broadband PSD plot figure handle
    end
    if PROCSPLD
        Param.figspld = figSplD; % decidecade SPL plot figure handle
    end

    if curChan == 1 % do not repeat info display for next channels
        if PROCPSDB || PROCSPLD
            fprintf(' FFT: %d points, %d%% overlap', Param.nf,FFTOVLAP);
            fprintf(', resolution = %.3f Hz\n', fVect(2));
            fprintf(' Snapshot contains %g FFT segments\n', Param.fn);
        end
        if HIPASSF > 0 || LOPASSF > 0
            fprintf(' Filter = [%.3f - %.3f] Hz\n',HIPASSF,LOPASSF);
        end
    end

    %%%%%%% Decidecade (DD) band setup: center f and limits
    if PROCSPLD
        DDmat = []; % initialize DD index matrix
        for idx = DDIDXL:DDIDXH
            % Decidecade center f are nominal ISO preferred, while
            % band limits are computed using power-of-ten formula.
            % Index 1 is for 1-Hz band, need to shift by 1 to get 10^0.
            DDfl = 10^( (2*(idx-1)-1) / 20 );
            DDfh = 10^( (2*(idx-1)+1) / 20 );
            fVidxl = find(fVect>=DDfl, 1, 'first'); % fvect low idx
            fVidxh = find(fVect<=DDfh, 1, 'last');  % fvect hi idx
            if DDfl > smpRate/2 || DDfh > smpRate/2
                exStat = 'Error: topmost decidecade too high.';
                disp(exStat)
                retVal = 1;
                msgbox(exStat,'Error','error','modal');
                return
            end
            if VERB > 0
                fprintf('DD(%.1f)=[%.3f-%.3f], [%d-%d]=[%.3f-%.3f]\n', ...
                    ISOF(idx),DDfl,DDfh, ...
                    fVidxl, fVidxh, fVect(fVidxl),fVect(fVidxh) );
            end           
            if fVect(2) > DDfh - DDfl % f resolution too coarse
                exStat='Error: increase FFT points to resolve DD band.';
                disp(exStat)
                retVal = 1;
                msgbox(exStat,'Error','error','modal');
                return
            end
            %%% DD matrix: ISO [f0, flow, fhi]; computed [flow, fhi];
            %   vector indexes [low, hi]
            DDmat = [DDmat; ISOF(idx) DDfl DDfh ...
                fVect(fVidxl) fVect(fVidxh) fVidxl fVidxh];
        end
        DDf0 = DDmat(:,1); % vector of ISO preferred center freq.s
        % interpolate calibration SF values in ISO center frequencies
        DDCalSF = interp1(fVect,CalSF,DDf0,'linear');
    else
        % pass void arguments to subroutine if DD not selected
        DDmat = []; DDCalSF = []; 
    end %%%%%%%%%% Decidecade band setup
    
    % initialize each results matrix
    outMatSplB = []; outMatPsdB = []; outMatSplD = []; outMatSplP = [];
    
    % reduce figure height for small screens
    scrSiz = get(0,'ScreenSize');
    if FIGSIZE(2)+FIGXYPOS(2) > scrSiz(4) % figure top goes out of screen
        FIGSIZE(2) = scrSiz(4)-FIGXYPOS(2);
    end

    % clear and reset position of all figures if requested
    clf(figWav);
    if PROCSPLB
        clf(figSplB)
    end
    if PROCSPLP
        clf(figSplP)
    end
    if PROCPSDB
        clf(figPsdB)
    end
    if PROCSPLD
        clf(figSplD)
    end
    if ~FREEFIG
        set(0,'CurrentFigure',figWav);
        figWav.OuterPosition=[FIGXYPOS FIGSIZE];
        if PROCSPLB
            set(0,'CurrentFigure',figSplB);
            figSplB.OuterPosition=[FIGXYPOS+FIGXYOFFS FIGSIZE];
        end
        if PROCPSDB
            set(0,'CurrentFigure',figPsdB);
            figPsdB.OuterPosition=[FIGXYPOS+2*FIGXYOFFS FIGSIZE];
        end
        if PROCSPLD
            set(0,'CurrentFigure',figSplD);              
            figSplD.OuterPosition=[FIGXYPOS+3*FIGXYOFFS FIGSIZE];
        end
        if PROCSPLP
            set(0,'CurrentFigure',figSplP);              
            figSplP.OuterPosition=[FIGXYPOS+4*FIGXYOFFS FIGSIZE];
        end
    end
    if DISPINFO % set info figure contents and properties
        infoTx1=sprintf(['WAV: %s\n    sr = %.1f kHz,    ' ...
            'len = %.1f s\n\n'], WavFileName, smpRate/1e3, fileDur);
        infoTx2=sprintf(['CAL: %s\n    %d pts, (%g-%g) Hz, ' ...
            'mean SF = %.1f Pa\n\n'], [CalFileName CalFileExt], ...
            calPts, calMat(1,1), calMat(calPts,1), meanCalSF );
        infoTx3=sprintf('FFT res. = %.2f Hz\n\n', fVect(2) );
        if nSnap > 1
            infoTx4=sprintf('%d snaps,    %.2f FFT segments each\n', ...
                nSnap, Param.fn);
        else
            infoTx4=sprintf('1 snap,    %.2f FFT segments\n', Param.fn);
        end
        figInfo = msgbox([infoTx1 infoTx2 infoTx3 infoTx4], ...
            'Info','replace'); % redraws on same box for each file
        figInfoPos = get(figInfo,'Position');
        set(figInfo,'Position',[FIGINFOXY figInfoPos(3) figInfoPos(4)]);
    end

    if nChans > 1 % stereo or multichannel
        fprintf(' Processing channel %d of %d:\n',curChan,nChans);
    end

    %%%%%%%%%%%%%%%%%%%%%%% snapshot loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for snapN = 1 : nSnap

        Param.sn = snapN;  % passed on to subroutines

        smpToRead = snapSmp(snapN,:); % 2-element vect [start,stop]
        fprintf(' Snap %d/%d t=[%.1f-%.1f] s=[%d-%d] ',...
            snapN,nSnap,smpToRead/smpRate,smpToRead);
        
        % read snapshot from file
        [AllChanY,smpRate] = audioread(WavFile,smpToRead);

        if nChans > 1
            OrigY = AllChanY(:,curChan); % select current channel
        else
            OrigY = AllChanY; % mono
        end

        Snap.Y = OrigY - mean(OrigY); % remove DC component
        snapStart = smpToRead(1)/smpRate; % update start time

        %%%%% plot snap waveform first anyway
        figure(figWav); % on top
        % set(0,'CurrentFigure',figWav);
        plot(Snap.T + snapStart, Snap.Y * meanCalSF)
        grid on
        xlabel(FigLab.XLABWAV)
        ylabel(FigLab.YLABWAV)
        if TITL
            title(sprintf(FigLab.TITLWAV,WavFileName,snapN,nSnap), ...
              'Interpreter','none') % otherwise underscores in filename
                                    % are interpreted as subscript
        else
            title('')
        end
        if HIPASSF > 0 || LOPASSF > 0
            % create annotation box, filter, replot filtered wave
            figure(FIGWAV); % on top
            pause(PLTPAUSE) % show unfiltered wave for a while
            ann = annotation('textbox',FILTBOXPOS, ...
                'string','Filtering ...');
            ann.HorizontalAlignment = 'center';
            ann.VerticalAlignment = 'middle';
            ann.BackgroundColor = [1 1 1];
            fprintf(' filt... '); % takes time
            pause(0.5) % need this to correctly show annotation
            if HIPASSF == 0
                if LOPASSF == 0
                    % no filter
                else
                    Snap.Y = lowpass(OrigY, LOPASSF, smpRate); 
                end
            else % highpass f is > 0
                if LOPASSF == 0
                    Snap.Y = highpass(OrigY, HIPASSF, smpRate); 
                else
                    Snap.Y = bandpass(OrigY, [HIPASSF LOPASSF], smpRate); 
                end
            end
            % default atten. 60dB, ripple 0.1dB, trans.width = 15%*f_cut
            clf; % needed to clear annotation box
            plot(Snap.T + snapStart, Snap.Y * meanCalSF)
            grid on % re-apply grid, labels and title as above
            xlabel(FigLab.XLABWAV)
            ylabel(FigLab.YLABWAV)
            if TITL
                title(sprintf([FigLab.TITLWAV ' filtered'], ...
                    WavFileName,snapN,nSnap), 'Interpreter','none')
            else
                title('')
            end
        end % filter
        pause(PLTPAUSE) 

        %%%%%%%%%%%% PSD subroutine call, processing in BB and/or DD

        if PROCPSDB || PROCSPLD || PROCSPLB || PROCSPLP
            [curPsdB,curSplD,curSplB,curSplP] = ...
                    TUNEsub41(Snap,FFTW,CalSF,DDCalSF,DDmat,Param);

            outMatPsdB = [outMatPsdB curPsdB];
            outMatSplD = [outMatSplD curSplD];
            outMatSplB = [outMatSplB curSplB]; % it is actually a vector
            outMatSplP = [outMatSplP curSplP]; %   <as above>
            outdbMatSplB   = 10*log10(outMatSplB)  - refdb;
            outdbMatSplP   = 20*log10(outMatSplP)  - refdb; % re 1 microPa
        end
        % pause already included at end of subroutine

        if PROCSPLB        % update broadband SPL plot for each snap
            figure(figSplB); % on top
            % set(0,'CurrentFigure',figSplB);
            % points are drawn along x at snapshot midpoints
            snapMidTimes = (snapTime(1:snapN,2)+snapTime(1:snapN,1)) / 2; 
            plot(snapMidTimes,outdbMatSplB,':ob','MarkerFaceColor','b')
            grid on
            xlabel(FigLab.XLABSPLB)
            ylabel(FigLab.YLABSPLB)
            ax = gca;  % get current axes handle
            ax.XScale = 'linear';
            xlim([STARTTIME stopTime]) % set time axis to full interval
            if YAXAUTO
                ylim auto % autorange: override Y limits
            else
                ylim([YAXMIN YAXMAX])
            end
            if TITL
                title(sprintf(FigLab.TITLSPLB, ...
                    WavFileName,STARTTIME,nSnap*snapLen,snapLen), ...
                    'Interpreter','none') 
            else
                title('')
            end
        end % BB SPL plot

        if PROCSPLP        % update zero-to-peak SPL plot for each snap
            figure(figSplP); % on top
            % set(0,'CurrentFigure',figSplB);
            % points are drawn along x at snapshot midpoints
            snapMidTimes = (snapTime(1:snapN,2)+snapTime(1:snapN,1)) / 2; 
            plot(snapMidTimes,outdbMatSplP,':ob','MarkerFaceColor','b')
            grid on
            xlabel(FigLab.XLABSPLP)
            ylabel(FigLab.YLABSPLP)
            ax = gca;  % get current axes handle
            ax.XScale = 'linear';
            xlim([STARTTIME stopTime]) % set time axis to full interval
            if YAXAUTO
                ylim auto % autorange: override Y limits
            else
                ylim([YAXMIN YAXMAX])
            end
            if TITL
                title(sprintf(FigLab.TITLSPLP, ...
                    WavFileName,STARTTIME,nSnap*snapLen,snapLen), ...
                    'Interpreter','none') 
            else
                title('')
            end
        end % peak SPL plot

        fprintf('\n'); % newline in command window when done

    end % snapshot loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% BB, DD averaging (done on linear values, not dB)
    switch avgType 
        case 1
            avgStr = 'Arithmetic mean';
            outAvgPsdB = mean(outMatPsdB,2);
            outAvgSplD = mean(outMatSplD,2);
        case 2
            avgStr = 'Median';
            outAvgPsdB = median(outMatPsdB,2);
            outAvgSplD = median(outMatSplD,2);
        otherwise % default: should not get here anyway
            avgStr = 'Arithmetic mean';
            outAvgPsdB = mean(outMatPsdB,2);
            outAvgSplD = mean(outMatSplD,2);
    end

    %%%%% percentiles for psd
    outPcHPsdB = prctile(outMatPsdB,PRCTILEH,2);
    outPcLPsdB = prctile(outMatPsdB,PRCTILEL,2);
    outPcHSplD = prctile(outMatSplD,PRCTILEH,2);
    outPcLSplD = prctile(outMatSplD,PRCTILEL,2);

    %%%%% convert results to dB re 1 microPa
    outdbMatPsdB   = 10*log10(outMatPsdB) - refdb; % full matrix 
    outdbAvgPsdB   = 10*log10(outAvgPsdB) - refdb; % average vector
    outdbPctHPsdB  = 10*log10(outPcHPsdB) - refdb;
    outdbPctLPsdB  = 10*log10(outPcLPsdB) - refdb;
    outdbMatSplD   = 10*log10(outMatSplD) - refdb; % full matrix 
    outdbAvgSplD   = 10*log10(outAvgSplD) - refdb; % average vector
    outdbPctHSplD  = 10*log10(outPcHSplD) - refdb;
    outdbPctLSplD  = 10*log10(outPcLSplD) - refdb;

    %%%%%%%%%%%%%%%%%%%%%% plot results according to processing type
    if PROCSPLB || PROCSPLP
        % no new plot: already updated on each snap
    end
    
    if PROCPSDB         % PSD
        set(0,'CurrentFigure',figPsdB);
        clf % do not redraw on last snapshot plot
        hold on
        if DRAWALL > 0
            for snapN = 1 : nSnap % draw all snapshots
               semilogx( fVect, outdbMatPsdB(:,snapN), 'c' )
            end
        end
        % return values are used to reference legend items
        psdL = semilogx( fVect, outdbPctLPsdB, 'g' ); % low percentile
        psdH = semilogx( fVect, outdbPctHPsdB, 'r' ); % hi percentile
        psdA = semilogx( fVect, outdbAvgPsdB, 'b' ); % avg: on top
        grid on
        hold off
    end
    
    if PROCSPLD         % decidecade SPL
        % line plot: selected average type
        % boxplot: median, fixed 25%-75% percentiles, min-max, outliers
        % transpose matrix, series must be column-wise
        figure(figSplD); % on top
        % set(0,'CurrentFigure',figSplD);
        if BOXPLT 
            if nSnap == 1 % no statistics for percentiles, use single plot
                fprintf('\n Single snapshot: not enough for boxplot\n');
                semilogx(DDf0,outdbAvgSplD,'ob', ...
                    'MarkerFaceColor','b','MarkerSize',8)
            else % plot both avg bullets and boxplot with percentiles
                plot(1:length(DDf0),outdbAvgSplD,'ob', ...
                    'MarkerFaceColor','b','MarkerSize',8)
                hold on
                boxplot(outdbMatSplD','Labels',num2str(DDf0))
                hold off
            end
        else % plot three bullets for low-hi percentiles and average
            ddL = semilogx(DDf0,outdbPctLSplD,'og');
            hold on
            ddH = semilogx(DDf0,outdbPctHSplD,'or');
            ddA = semilogx(DDf0,outdbAvgSplD,'ob','MarkerFaceColor','b');
            hold off
        end
    end

    %%%%% overlay Knudsen curves for N sea states - only for psd
    if PROCPSDB && KNDSDRAW
        figure(figPsdB); % on top
        % set(0,'CurrentFigure',figPsdB);
        hold on
        seaStates = SSLABELS; 
        kndsIdx = KNDSSIDX + 1; % SS range from 0, index from 1
        for idx = 1:length(kndsIdx)
            knudsNL1K = KNDLEVS(kndsIdx(idx));
            knudsNL = knudsNL1K + KNDDEC*log10(fVect/1e3); 
            txtYpos = knudsNL1K + KNDDEC*log10(fVect(1)/1e3); 
            for idx2 = 1:length(fVect)
                if fVect(idx2) < KNDCUT
                    knudsNL(idx2) = knudsNL1K;
                    txtYpos = knudsNL1K;
                end
            end
            semilogx( fVect, knudsNL, '--m' )
            text(TXTXPOS,txtYpos+TXTYOFFS,seaStates(kndsIdx(idx),:))
        end
        hold off
    end

    %%%%%%%%% apply type-specific title, scale, labels, opt. legend.
    % Note: interpreter set to 'none' in title otherwise underscores
    %       in filename are interpreted as subscript.

    if PROCSPLB || PROCSPLP
        % plot already complete
    end

    if PROCPSDB
        figure(figPsdB); % on top
        % set(0,'CurrentFigure',figPsdB); % update without focus
        grid on
        xlabel(FigLab.XLABPSDB)
        ylabel(FigLab.YLABPSDB)
        ax = gca;  % get current axes handle
        ax.XScale = 'log';
        if XAXAUTO
            xlim auto % autorange: override X limits
        else
            xlim([XAXMIN XAXMAX])
        end
        if YAXAUTO
            ylim auto % autorange: override Y limits
        else
            ylim([YAXMIN YAXMAX])
        end
        if TITL
            title(sprintf(FigLab.TITLPSDB,...
                WavFileName,STARTTIME,nSnap*snapLen,snapLen,NFFT), ...
                'Interpreter','none')
        else
            title('')
        end
        if LEGND
            legend( [psdH psdA psdL], ...
                FigLab.LGNDPSDB(2), avgStr, FigLab.LGNDPSDB(1) )
        end
    end

    if PROCSPLD
        figure(figSplD); % on top
        % set(0,'CurrentFigure',figSplD); % update without focus
        grid on
        xlabel(FigLab.XLABSPLD)
        ylabel(FigLab.YLABSPLD)
        ax = gca;  % get current axes handle
        if BOXPLT
            if nSnap == 1
                ax.XScale = 'log'; % no stats for percentiles
            else
                ax.XScale = 'linear'; % boxplot decidecade series
            end
        else
            ax.XScale = 'log'; % decidecade line plot
        end
        ylim([YAXMIN YAXMAX]) % only resize Y axis, X set by DD bands
        if YAXAUTO
            ylim auto % autorange: override Y limits
        end
        if TITL
            title(sprintf(FigLab.TITLSPLD,...
                WavFileName,STARTTIME,nSnap*snapLen,snapLen,NFFT), ...
                'Interpreter','none') 
        else
            title('')
        end
        if LEGND
            if BOXPLT && nSnap > 1 % No boxplot output for single snap.
                % Use trick to split legend in multiple rows to
                % align linetype mark with string of average type;
                % rows must be equal length, pad avgStr with spaces
                spPad = length(FigLab.LGNDBOX)-strlength(avgStr);
                legend([ FigLab.LGNDBOX(1,:); ...
                         FigLab.LGNDBOX(2,:); ...
                         [avgStr char(32*ones(1,spPad))]; ...
                         FigLab.LGNDBOX(3,:); ...
                         FigLab.LGNDBOX(4,:) ]) 
            end
            if ~BOXPLT % decidecade line plot, also for single snap
                legend( [ddH ddA ddL], ...
                FigLab.LGNDSPLD(2), avgStr, FigLab.LGNDSPLD(1) )

            end
        end
    end

    % Add date_time in output filenames (jpeg&text) to make them unique
    dtStr = string(datetime('now'), '_yyMMdd_HHmm');
    
    %%%%% Save figures as *.jpg file with same size as on screen
    fullPath = [outPath FileNameNoExt]; % File path and name without ext.
    if SAVEJPG
        disp([' Saving figures to ' outPath ' :'])
        if PROCSPLB
            set(FIGSPLB,'PaperPositionMode','auto')
            jpegFile = strcat(fullPath, dtStr, '_bbspl.jpg');
            print(['-f' num2str(FIGSPLB)],FIGDRV,FIGRES,jpegFile)
            [~,fName,fExt] = fileparts(jpegFile);
            disp(' - ' + fName + fExt)
        end
        if PROCSPLP
            set(FIGSPLP,'PaperPositionMode','auto')
            jpegFile = strcat(fullPath, dtStr, '_pkspl.jpg');
            print(['-f' num2str(FIGSPLP)],FIGDRV,FIGRES,jpegFile)
            [~,fName,fExt] = fileparts(jpegFile);
            disp(' - ' + fName + fExt)
        end
        if PROCPSDB
            set(FIGPSDB,'PaperPositionMode','auto')
            jpegFile = strcat(fullPath, dtStr, '_bbpsd.jpg');
            print(['-f' num2str(FIGPSDB)],FIGDRV,FIGRES,jpegFile)
            [~,fName,fExt] = fileparts(jpegFile);
            disp(' - ' + fName + fExt)
        end
        if PROCSPLD
            set(FIGSPLD,'PaperPositionMode','auto')
            if BOXPLT
                jpegFile = strcat(fullPath, dtStr, '_ddbox.jpg');
                print(['-f' num2str(FIGSPLD)],FIGDRV,FIGRES,jpegFile)
                [~,fName,fExt] = fileparts(jpegFile);
                disp(' - ' + fName + fExt)
            else
                jpegFile = strcat(fullPath, dtStr, '_ddspl.jpg');
                print(['-f' num2str(FIGSPLD)],FIGDRV,FIGRES,jpegFile)
                [~,fName,fExt] = fileparts(jpegFile);
                disp(' - ' + fName + fExt)
            end
        end
    end

    %%%%%% Save results as text files.
    % Transpose results vectors to build output matrix,
    % use rounding to reduce n. of decimals from default = 13.
    % Note: function dlmwrite() not recommended after R2019,
    %       use function writematrix() instead
    if SAVETXT
        disp([' Saving text files to ' outPath ' :'])
        if PROCSPLB
            outTxtMat = [snapMidTimes,outdbMatSplB'];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_bbspl.txt');
            outTxtId = fopen(outTxtFile,'wt'); 
            fprintf(outTxtId,HDRSPLB); % write header
            fclose (outTxtId);
            if versYear < 2019 % must use old function (not recommended)
                dlmwrite(outTxtFile,outTxtMat,'-append','delimiter','\t') 
            else % use new function introduced in 2019
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab',...
                                                 'WriteMode','append')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(' - ' + fName + fExt)
        end

        if PROCSPLP
            outTxtMat = [snapMidTimes,outdbMatSplP'];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_pkspl.txt');
            outTxtId = fopen(outTxtFile,'wt'); 
            fprintf(outTxtId,HDRSPLP); % write header
            fclose (outTxtId);
            if versYear < 2019 % must use old function (not recommended)
                dlmwrite(outTxtFile,outTxtMat,'-append','delimiter','\t') 
            else % use new function introduced in 2019
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab',...
                                                 'WriteMode','append')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(' - ' + fName + fExt)
        end

        if PROCPSDB
            outTxtMat = [fVect',outdbAvgPsdB,outdbPctLPsdB,outdbPctHPsdB];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_bbpsd.txt');
            outTxtId = fopen(outTxtFile,'wt'); 
            fprintf(outTxtId,HDRPSDB); % write header
            fclose (outTxtId);
            if versYear < 2019 % must use old function (not recommended)
                dlmwrite(outTxtFile,outTxtMat,'-append','delimiter','\t') 
            else % use new function introduced in 2019
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab',...
                                                 'WriteMode','append')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(' - ' + fName + fExt)
        end

        if PROCSPLD
            outTxtMat = [DDf0,outdbAvgSplD,outdbPctLSplD,outdbPctHSplD];
            outTxtMat = round(outTxtMat*10^TXTDECS)/10^TXTDECS;
            outTxtFile = strcat(fullPath, dtStr, '_ddspl.txt');
            outTxtId = fopen(outTxtFile,'wt'); 
            fprintf(outTxtId,HDRSPLD); % write header
            fclose (outTxtId);
            if versYear < 2019 % must use old function (not recommended)
                dlmwrite(outTxtFile,outTxtMat,'-append','delimiter','\t') 
            else % use new function introduced in 2019
                writematrix(outTxtMat,outTxtFile,'Delimiter','tab',...
                                                 'WriteMode','append')
            end
            [~,fName,fExt] = fileparts(outTxtFile);
            disp(' - ' + fName + fExt)
        end
    end
    pause(PLTPAUSE * 4) % processing end, longer pause

    end % channels loop

end %%%%%%%%%% datafile loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pass existing figure handles to GUI so they can be closed on exit
if exist('figWav','var') > 0
    figHandles = figWav; % Should exist unless an error occurs
end
if exist('figInfo','var') > 0
    figHandles = [figHandles figInfo];
end
if exist('figSplB','var') > 0
    figHandles = [figHandles figSplB];
end
if exist('figSplP','var') > 0
    figHandles = [figHandles figSplP];
end
if exist('figPsdB','var') > 0
    figHandles = [figHandles figPsdB];
end
if exist('figSplD','var') > 0
    figHandles = [figHandles figSplD];
end

timeOfEnd = datetime('now');
elapsTime = sprintf(' in %s', timeOfEnd - timeOfStart);
switch validWFiles
    case 1
        exStat = 'Normal end - Processed 1 file';
    otherwise
        exStat = sprintf('Normal end - Processed %d files',validWFiles);
end
exStat = [exStat elapsTime]; % Add elapsed processing time

if retVal > 0
    exStat = [exStat ' - ERROR: see logfile'];
end

if LOGTOFILE
    if LOGSETUP % Display setup so it can also be logged to file
        fprintf('\nSetup:\n');
        disp(Setup)
    end
    fprintf('\nDone: %s\n\n',datetime('now'))
    diary off % Stop logging to file before exit, if enabled
else
    disp('Done.'); % Normal end
end

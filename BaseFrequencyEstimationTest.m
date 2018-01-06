%This programs compares how well a statistical Zero Crossing Frequency Estimation does ###########################################################################
% in comparison to one based on a FFT#############################################################################################################################


% Naming Conventions: ============================================================================================================================================
% Small letters are used before other abbreviatons and/or words starting with large letters
% Key variable names like 'fs' for 'sampling frequency' usually don't follow this pattern


% a : average
% b : beginning
% c : current
% d : defective
% diff : difference
% e : end
% f : first
% i : iterator
% l : last
% max : maximum
% min : minimum
% m : median
% n : number
% o : observe (after this ZC, the next ZC gets repeated, regardless of polarity)
% p : position
% r : relevant
% s : size
% t : time
% v : verified


% Cyc : Cycle
% Pol : Polarity
% Seg : Segment
% Win : Window
% ZC : Zero crossing

% Main Function ==================================================================================================================================================
function[] = SampleTest()
% Initializations ________________________________________________________________________________________________________________________________________________
  format long
  global fs;
  global DataOutputDir

  FigOutputDir = 'WaveformOutput\';           % Directories
  DataOutputDir = 'DataOutput\';
  mkdir (FigOutputDir);
  mkdir (DataOutputDir);

  PlotSignal = 1;           % set to 1 to print out signal waveform

  sWindow=1024;
  nFFTPoints = 1024;           % FFT-window
  nWindows = 100;
  % System conditions
  fs=44100; %Sampling frequency
  fc = 20000;           % Filter cutoff frequency
  f = 500; %Frequency of signal in Hz
  ExactFreq = 500;
  t = 0:1/fs:(((sWindow*nWindows) - 1)/fs); %time index
  Saw=sawtooth(2*pi*f*t);

  % Signal generation (Adding and Filtering)
  [b,a] = butter(2,fc/(fs/2));          % butterworth filter coefficients
  bSignal = filter(b,a,Saw);

  SNRdB = [90, 80, 70, 50, 45, 40, 35, 30, 25, 20, 18, 16, 14, 12, 10, 8, 6, 4, 2, 0];
  nSNR = Size(SNRdB);

  [WinSumOffsetAll, FFTWinSumOffsetAll] = InitString();           % OutputStrings
  [WinSumOffsetAllRedux, FFTWinSumOffsetAllRedux] = InitReduxString();

  [Counter, NaNCounter, nError, PrintError] = AssignVal2Var(0);

  sbSignal = Size(bSignal);
  nFullWindow = fix(sbSignal / sWindow);
  csWindow = sWindow;

% Main loop ______________________________________________________________________________________________________________________________________________________
  for cSNR=1:nSNR           % SNR-loop
    [Offset, WinSumOffset, FFTOffset, FFTWinSumOffset] = InitString();
    [OffsetRedux, WinSumOffsetRedux, FFTOffsetRedux, FFTWinSumOffsetRedux] = InitReduxString();
    [FFTCritNum, CritNum] = AssignVal2Var(0);

    bSignal = awgn(bSignal,SNRdB(cSNR));
    [PrintError] = PrintOutput(strcat(num2str(SNRdB(cSNR)),'Signal.data'), num2str(bSignal), PrintError);

    for cWindow=1:(nFullWindow)           % Window-loop
      ClearVars = {'pPosZC', 'pNegZC', 'nZC', 'FirstZC', 'STATE', 'prNegZC', 'prPosZC', 'rNegCyc', 'rPosCyc', 'drNegZC', 'drPosZC', 'odrPosZC', 'odrNegZC', 'oodrNegZC', 'oodrPosZC', 'FFT', 'FFT_temp', 'FFTFigXAxis', 'MedianArray'};
      clear(ClearVars{:});

      [bReSizer, eReSizer, csWindow] = WinSize(cWindow, sWindow, sbSignal);
      cSignal = bSignal(bReSizer:eReSizer);

      [pPosZC,pNegZC,FirstZC,lZC,nZC] = ZCD(cSignal, sWindow);          % lZC and nZC are unused as of now

      nPosZC = Size(pPosZC);
      nNegZC = Size(pNegZC);

% Plot and output signal waveform ________________________________________________________________________________________________________________________________
    if (PlotSignal == 1)            % default: 1 --> plot signal waveform
      if (cWindow == 50)            % middle window
        for iCounter=1:350
          Counter(iCounter) = iCounter;
        end

        gSignal = zeros(1,sWindow);
        sgSignal = size(gSignal);
        gSignal = bSignal(bReSizer:bReSizer + 349);

        Normalizer = max(abs(gSignal));
        gSignal = gSignal / (Normalizer + 0.0000001);

        FigOutputName = 'SNR_TestSaw';
        FigOutput = strcat(FigOutputDir, num2str(SNRdB(cSNR)), 'dB', FigOutputName);
        f = figure('visible', 'off');
        plot(Counter,gSignal);
        hold on;
        stem(Counter,gSignal);
        hold off;
        xlabel('Sample Number');
        ylabel('Sample Value');
        set(gca,'ytick',-1.2:0.2:1.2)

        saveas(f, FigOutput, 'png')
        saveas(f, FigOutput, 'fig')
        % matlab2tikz(strcat(FigOutput, '.tex'), 'showInfo', false);
      end
    end

% FFT basefreq estimation ________________________________________________________________________________________________________________________________________
      DFTBinFreq = fs / nFFTPoints;
      HalfFFT = round(nFFTPoints / 2);
      FFT = fft(cSignal);
      for j=1:nFFTPoints
        FFT(j) = FFT(j) / nFFTPoints;
      end
      FFT_temp = abs(FFT(1:round(nFFTPoints / 2)));

      [NOTUSED, FFTFreqIndex] = max(FFT_temp);
      FFTFreq(cWindow) = DFTBinFreq * (FFTFreqIndex - 1);

      FFTAbsoluteOffset(cWindow) = FFTFreq(cWindow) - ExactFreq;            % print data out
      FFTRelativeOffset(cWindow) = abs((FFTAbsoluteOffset(cWindow) / ExactFreq) * 100);
      if (FFTRelativeOffset(cWindow) > 5.0)
        FFTOffset = strcat(FFTOffset, '!!,,','Relative Error:,',num2str(FFTRelativeOffset(cWindow)),' percent, Absolute Error:,',num2str(FFTAbsoluteOffset(cWindow)),',Calculated Frequency:,',num2str(FFTFreq(cWindow)),'\n\n\n');
        FFTOffsetRedux = strcat(FFTOffsetRedux, num2str(FFTRelativeOffset(cWindow)),',,',num2str(FFTAbsoluteOffset(cWindow)),',,',num2str(FFTFreq(cWindow)),'\n\n');
        FFTCritNum = FFTCritNum + 1;
      elseif (FFTRelativeOffset(cWindow) <= 5.0)
        FFTOffset = strcat(FFTOffset,'Relative Error:,',num2str(FFTRelativeOffset(cWindow)),' percent, Absolute Error:,',num2str(FFTAbsoluteOffset(cWindow)),'Calculated Frequency:,',num2str(FFTFreq(cWindow)),'\n\n\n');
        FFTOffsetRedux = strcat(FFTOffsetRedux, num2str(FFTRelativeOffset(cWindow)),',,',num2str(FFTAbsoluteOffset(cWindow)),',,',num2str(FFTFreq(cWindow)),'\n\n');
      elseif (FFTRelativeOffset(cWindow) ~= FFTRelativeOffset(cWindow))
        disp(FFTRelativeOffset(cWindow))
        error('Error!\nFFTRelativeOffset is NaN:');
      else
        disp(FFTRelativeOffset(cWindow))
        error('Critical Error!\nInternal Logic Error #2');
      end

% ZC basefreq estimation _________________________________________________________________________________________________________________________________________
      lPosZC = nPosZC - 1;
      lNegZC = nNegZC - 1;

      if nPosZC == nNegZC           % necessary to establish order and limit of value-checks
        if FirstZC == 0
          [NegCyc, PosCyc, sNegSeg, sPosSeg] = sZCSearch((nPosZC-1), (nPosZC-1), pNegZC, pPosZC);
        else
          [PosCyc, NegCyc, sPosSeg, sNegSeg] = sZCSearch((nPosZC-1), (nPosZC-1), pPosZC, pNegZC);
        end
      else
        if FirstZC == 0            % more nNegZC
          [NegCyc, PosCyc, sNegSeg, sPosSeg] = sZCSearch((nNegZC - 1), (nPosZC - 1), pNegZC, pPosZC);
          else
          [PosCyc, NegCyc, sPosSeg, sNegSeg] = sZCSearch((nPosZC - 1), (nNegZC - 1), pPosZC, pNegZC);
        end
      end

      nPosCyc = Size(PosCyc);
      nNegCyc = Size(NegCyc);

      if PosCyc(1) > 0
        [smaxPosZC,smaxPosIndex] = max(PosCyc);
        if (FirstZC == 0)
          smaxPosSeg = pNegZC(smaxPosIndex + 1) - pPosZC(smaxPosIndex);
        else
          smaxPosSeg = pNegZC(smaxPosIndex) - pPosZC(smaxPosIndex);
        end
      else
        error('PosCyc is 0 for some reason')
      end

      if NegCyc(1) > 0
        [smaxNegZC,smaxNegIndex] = max(NegCyc);
        if (FirstZC == 0)
          smaxNegSeg = pPosZC(smaxNegIndex) - pNegZC(smaxNegIndex);
        else
          smaxNegSeg = pPosZC(smaxNegIndex + 1) - pNegZC(smaxNegIndex);
        end
      else
        error('NegCyc is 0 for some reason')
      end

      if ((smaxPosZC < 1) || (smaxNegZC < 1))
        error('error calculating smaxZC');
      end

      [rPosCyc, prPosZC, nrPosCyc, lrPosZC] = FindRelevantZCs(PosCyc, smaxPosZC, pPosZC, lPosZC, sPosSeg, smaxPosSeg);
      [rNegCyc, prNegZC, nrNegCyc, lrNegZC] = FindRelevantZCs(NegCyc, smaxNegZC, pNegZC, lNegZC, sNegSeg, smaxNegSeg);

      rSegSTATE = 0;
      STATE = 1;            % Used to prevent infinite loop

% repeat if no rdPolZC are found ---------------------------------------------------------------------------------------------------------------------------------
      if (isnan(nrPosCyc))           % repeat search for rPosZC
        i = 2;
        while (STATE ~= 0)
          if (i>nPosCyc)
            rSegSTATE = rSegSTATE + 1;
            STATE = 0;
          end
          [smaxPosZC, smaxPosIndex] = kLargestElement(PosCyc, i);
          if (FirstZC == 0)
            smaxPosSeg = pNegZC(smaxPosIndex + 1) - pPosZC(smaxPosIndex);
          else
            smaxPosSeg = pNegZC(smaxPosIndex) - pPosZC(smaxPosIndex);
          end
          [rPosCyc, prPosZC, nrPosCyc, lrPosZC] = FindRelevantZCs(PosCyc, smaxPosZC, pPosZC, lPosZC, sPosSeg, smaxPosSeg);
          if (isnan(nrPosCyc))
            i = i + 1;
          else
            STATE = 0;
          end
        end
      end

      if (isnan(nrNegCyc))
        STATE = 1;
        i = 2;
        while (STATE ~= 0)
          if (i>nNegCyc)
            rSegSTATE = rSegSTATE + 2;
            STATE = 0;
          end
          [smaxNegZC, smaxNegIndex] = kLargestElement(NegCyc, i);
          if (FirstZC == 0)
            smaxNegSeg = pPosZC(smaxNegIndex) - pNegZC(smaxNegIndex);
          else
            smaxNegSeg = pPosZC(smaxNegIndex + 1) - pNegZC(smaxNegIndex);
          end
          [rNegCyc, prNegZC, nrNegCyc, lrNegZC] = FindRelevantZCs(NegCyc, smaxNegZC, pNegZC, lNegZC, sNegSeg, smaxNegSeg);
          if (isnan(nrNegCyc))
            i = i + 1;
          else
            STATE = 0;
          end
        end
      end

      RelativeOffset = NaN;
      if (rSegSTATE == 3)           % no rdPolZCs found at all
        AbsoluteOffset(cWindow) = NaN;
        RelativeOffset(cWindow) = NaN;
        Offset = strcat(Offset, '!! No usable rZCs found at all\n\n\n');
        OffsetRedux = strcat(OffsetRedux, 'No usable rZCs found at all\n\n\n');
        CritNum = CritNum + 1;
      else
        if (rSegSTATE == 1)           % only multiple rdNegZCs found
          [N_NegCyc, aN_NegCyc, nN_NegCyc] = FindPolCyc(prNegZC, nrNegCyc);
          [Freq(cWindow), nPeriod] = SinglePolStatfilter(N_NegCyc);
          nFFTPoints = round(nPeriod);
          DFTBinFreq = fs / nFFTPoints;

          AbsoluteOffset(cWindow) = Freq(cWindow) - ExactFreq;
          RelativeOffset(cWindow) = abs((AbsoluteOffset(cWindow) / ExactFreq) * 100);
          Offset = strcat(Offset, 'No rPosZCs ,');
        elseif (rSegSTATE == 2)           % only multiple rdPosZCs found
          [N_PosCyc, aN_PosCyc, nN_PosCyc] = FindPolCyc(prPosZC, nrPosCyc);
          [Freq(cWindow), nPeriod] = SinglePolStatfilter(N_PosCyc);
          nFFTPoints = round(nPeriod);
          DFTBinFreq = fs / nFFTPoints;

          AbsoluteOffset(cWindow) = Freq(cWindow) - ExactFreq;
          RelativeOffset(cWindow) = abs((AbsoluteOffset(cWindow) / ExactFreq) * 100);
          Offset = strcat(Offset, 'No rNegZCs ,');
        else            % both rdPolZCs found
          if((nrNegCyc < 2) || (nrPosCyc < 2))
            error('this should never occur, but too little nrsZCs')
          end
          [frZC,pfrZC,lrZC,plrZC] = FindFirstAndLastRelevantZCs(rPosCyc(1), rNegCyc(1), lrPosZC, lrNegZC);
          [N_PosCyc, aN_PosCyc, nN_PosCyc] = FindPolCyc(prPosZC, nrPosCyc);
          [N_NegCyc, aN_NegCyc, nN_NegCyc] = FindPolCyc(prNegZC, nrNegCyc);
          [nPeriod, Freq(cWindow)] = StatFilter (N_PosCyc, N_NegCyc);
        end

% Calculate error and print out (window) _________________________________________________________________________________________________________________________
        AbsoluteOffset(cWindow) = Freq(cWindow) - ExactFreq;            % print data out
        RelativeOffset(cWindow) = abs((AbsoluteOffset(cWindow) / ExactFreq) * 100);
        if (RelativeOffset(cWindow) > 5.0)
          Offset = strcat(Offset, '!!,,','Relative Error:,',num2str(RelativeOffset(cWindow)),' percent, Absolute Error:,',num2str(AbsoluteOffset(cWindow)),',Calculated Frequency:,',num2str(Freq(cWindow)),'\n\n\n');
          OffsetRedux = strcat(OffsetRedux, num2str(RelativeOffset(cWindow)),',,',num2str(AbsoluteOffset(cWindow)),',,',num2str(Freq(cWindow)),'\n\n');
          CritNum = CritNum + 1;
        elseif (RelativeOffset(cWindow) <= 5.0)
          Offset = strcat(Offset,'Relative Error:,',num2str(RelativeOffset(cWindow)),' percent, Absolute Error:,',num2str(AbsoluteOffset(cWindow)),'Calculated Frequency:,',num2str(Freq(cWindow)),'\n\n\n');
          OffsetRedux = strcat(OffsetRedux, num2str(RelativeOffset(cWindow)),',,',num2str(AbsoluteOffset(cWindow)),',,',num2str(Freq(cWindow)),'\n\n');
        elseif (RelativeOffset(cWindow) ~= RelativeOffset(cWindow))
          disp(RelativeOffset(cWindow))
          error('Error!\nRelativeOffset is NaN:');
        else
          disp(RelativeOffset(cWindow))
          error('Critical Error!\nInternal Logic Error #2');
        end
      end
    end

    [WinSumAbsoluteOffset, WinSumRelativeOffset, WinSumFreq, nSkip] = AssignVal2Var(0);
    WindowSkip = isnan(RelativeOffset);           % Prevent hard calculation errors from 'contaminating' mean-calculation.
    for cWindow=1:1
      if (WindowSkip(cWindow) == 0)
        WinSumFreq = WinSumFreq + Freq(cWindow);
      elseif (WindowSkip(cWindow) == 1)
        nSkip = nSkip + 1;
      else
        disp('Critical internal logic error, "isnan" function is not working!\n')
      end
    end

% Calculate sum of windows _______________________________________________________________________________________________________________________________________
% ZCFE part -------------------------------------------------------------------------------------------------------------------------------------------------------
    rWinSumWindows = nFullWindow - nSkip;
    WinSumFreq = mean(Freq); 
    WinSumAbsoluteOffset = WinSumFreq - ExactFreq;
    WinSumRelativeOffset = abs((WinSumAbsoluteOffset / ExactFreq) * 100);
    WinSumOffset = strcat('__AVERAGED__,,', 'Relative Error:, ',num2str(WinSumRelativeOffset),', Absolute Error:, ', num2str(WinSumAbsoluteOffset), ', True Frequency:, ',num2str(ExactFreq),', Calculated Frequency:, ', num2str(WinSumFreq), 'Times with more than 5 percent error:,', num2str(CritNum));
    WinSumOffsetRedux = strcat('AVERAGED:\n', num2str(WinSumRelativeOffset), ',,', num2str(WinSumAbsoluteOffset), ',,', num2str(ExactFreq), ',,', num2str(WinSumFreq), ',,', num2str(CritNum),'\n\n');

    if WinSumRelativeOffset >= 5
      WinSumOffset = strcat('!!,,' ,WinSumOffset);
    end
    WinSumOffset = strcat(WinSumOffset,'\n');
    WinSumOffsetAll = strcat(WinSumOffsetAll, num2str(SNRdB(cSNR)), 'dB:,', WinSumOffset);
    WinSumOffsetAllRedux = strcat(WinSumOffsetAllRedux, num2str(SNRdB(cSNR)), 'dB:,', WinSumOffsetRedux);

    [PrintError] = PrintOutput(strcat(num2str(SNRdB(cSNR)) ,'dB_ZCD.csv'), strcat(OffsetRedux, WinSumOffsetRedux), PrintError);
    [PrintError] = PrintOutput(strcat(num2str(SNRdB(cSNR)) ,'dB_ZCDVerbose.csv'), strcat(Offset, WinSumOffset), PrintError);

% FFT part --------------------------------------------------------------------------------------------------------------------------------------------------------
    [FFTWinSumAbsoluteOffset, FFTWinSumRelativeOffset, FFTWinSumFreq] = AssignVal2Var(0);
    FFTWinSumFreq = mean(FFTFreq);
    FFTWinSumAbsoluteOffset = FFTWinSumFreq - ExactFreq;
    FFTWinSumRelativeOffset = abs((FFTWinSumAbsoluteOffset / ExactFreq) * 100);
    FFTWinSumOffset = strcat('__AVERAGED__,,', 'Relative Error:, ',num2str(FFTWinSumRelativeOffset),', Absolute Error:, ', num2str(FFTWinSumAbsoluteOffset), 'Calculated Frequency:, ', num2str(FFTWinSumFreq), 'Times with more than 5 percent error:,', num2str(FFTCritNum));
    FFTWinSumOffsetRedux = strcat('AVERAGED:\n', num2str(FFTWinSumRelativeOffset), ',,', num2str(FFTWinSumAbsoluteOffset), ',,', num2str(FFTWinSumFreq), ',,', num2str(FFTCritNum),'\n\n');

    if FFTWinSumRelativeOffset >= 5
      FFTWinSumOffset = strcat('!!,,' ,FFTWinSumOffset);
    end
    FFTWinSumOffset = strcat(FFTWinSumOffset,'\n');
    FFTWinSumOffsetAll = strcat(FFTWinSumOffsetAll, num2str(SNRdB(cSNR)), 'dB:,', FFTWinSumOffset);
    FFTWinSumOffsetAllRedux = strcat(FFTWinSumOffsetAllRedux, num2str(SNRdB(cSNR)), 'dB:,', FFTWinSumOffsetRedux);

    [PrintError] = PrintOutput(strcat(num2str(SNRdB(cSNR)) ,'dB_FFTVerbose.csv'), strcat(FFTOffset, FFTWinSumOffset), PrintError);
    [PrintError] = PrintOutput(strcat(num2str(SNRdB(cSNR)) ,'dB_FFT.csv'), strcat(FFTOffsetRedux, FFTWinSumOffsetRedux), PrintError);


    ClearVars = {'AbsoluteOffset', 'RelativeOffset', 'WinSumAbsoluteOffset', 'WinSumRelativeOffset', 'Offset', 'WinSumOffset', 'FFTAbsoluteOffset', 'FFTRelativeOffset', 'FFTWinSumAbsoluteOffset', 'FFTWinSumRelativeOffset', 'FFTOffset', 'FFTWinSumOffset'};
    clear(ClearVars{:});
  end

% Print summaries to file _________________________________________________________________________________________________________________________________________
  [PrintError] = PrintOutput('FFTSummaryVerbose.csv', FFTWinSumOffsetAll, PrintError);
  [PrintError] = PrintOutput('FFTSummary.csv', FFTWinSumOffsetAllRedux, PrintError);
  [PrintError] = PrintOutput('ZCDSummaryVerbose.csv', WinSumOffsetAll, PrintError);
  [PrintError] = PrintOutput('ZCDSummary.csv', WinSumOffsetAllRedux, PrintError);

  if (PrintError ~= 0)
    disp(strcat('There have been ', num2str(PrintError), ' Errors writing data to file.'));
  else
    disp('There were no file writing errors.');
  end
end








% Functions ======================================================================================================================================================
function [Var1, Var2, Var3, Var4] = InitReduxString(); % Multiple Strings use this header
  Var1 = 'Relative error,,Absolute error,,Calculated Frequency,,More than 5 percent offset\n\n\n';
  Var2 = 'Relative error,,Absolute error,,Calculated Frequency,,More than 5 percent offset\n\n\n';
  Var3 = 'Relative error,,Absolute error,,Calculated Frequency,,More than 5 percent offset\n\n\n';
  Var4 = 'Relative error,,Absolute error,,Calculated Frequency,,More than 5 percent offset\n\n\n';
end

function [CycFirstZC, CycSecondZC, sFirstSeg, sSecondSeg] = sZCSearch(fLimit, lLimit, FirstZC, SecondZC) % Calculate distances from
  for j=1:fLimit       % Find all First ZC halfCycs
    CycFirstZC(j) = FirstZC(j+1) - FirstZC(j);
    sFirstSeg(j) = SecondZC(j) - FirstZC(j);
  end
  if lLimit >= 1
    for j=1:lLimit       % Find all Second Polarity halfCycs
      CycSecondZC(j) = SecondZC(j + 1) - SecondZC(j);
      sSecondSeg(j) = FirstZC(j + 1) - SecondZC(j);
    end
  else
    CycSecondZC(j) = 0;
  end
end

function [Cyc, pCyc] = sCycSearchSingle(Limit, ZC) % in case not enough nZC are found for two Pols
  for j=1:Limit 
    Cyc(j) = ZC(j+1) - ZC(j);
    pCyc(j) = ZC(j);
  end
end


function [FreqCyc] = SampleToFrequency(SumCyc, nCyc)
  global fs;
  sCyc = (SumCyc / nCyc) * (1.0 / fs);
  FreqCyc = 1 / sCyc;
end

function [pPosZC,pNegZC,FirstZC,lZC,nZC] = ZCD(cSignal, nPoints)
  nPoints;
  pNegZC(1) = 0.0;
  pPosZC(1) = 0.0;
  l=1;
  m=1;
  nZC = 0;
  if (cSignal(2) <= 0.0)
    STATE = 0;
    FirstZC = 1;
  else
    STATE = 1;
    FirstZC = 0;
  end
% find zero crossings _____________________________________________________________
  for j=2:nPoints
    if (cSignal(j) >= 0)        % positive Polarity
      if (STATE == 0)           % Polarity change?
        pPosZC(l) = j;
        nZC = nZC + 1;
        l = l + 1;
        lZC = 1;
        STATE = 1;
      end
    else                        % negative Polarity
      if (STATE == 1)           % Polarity change?
        pNegZC(m) = j;
        nZC = nZC + 1;
        m = m + 1;
        lZC = 0;
        STATE = 0;
      end
    end
  end
end

function [frZC,pfrZC,lrZC,plrZC] = FindFirstAndLastRelevantZCs(frPosCyc, frNegCyc, lrPosZC, lrNegZC)
  if (lrPosZC > lrNegZC)            % determine Polarity and position of last relevant ZCe
    lrZC = 1;
    plrZC = lrPosZC;
  else
    lrZC = 0;
    plrZC = lrNegZC;
  end

  if (frPosCyc < frNegCyc)            % case: nPos == nNeg
    frZC = 1;
    pfrZC = frPosCyc;
  else
    frZC = 0;
    pfrZC = frNegCyc;
  end
end


function [rsCyc, prCyc, nrsCyc, lrZC] = FindRelevantZCs(sCyc, maxsCyc, pZC, lZC, sSeg ,smaxSeg);
  m = 1;
  if maxsCyc > 0
    for j=1:lZC           % Find all windows of relevant length (maxlength / 2)
      if ((sCyc(j) >= (maxsCyc * 0.75)) && ((sSeg(j) >= (smaxSeg * 0.75)) && (sSeg(j) <= (smaxSeg * 1.25))))
        rsCyc(m) = sCyc(j);
        prCyc(m) = pZC(j);
        m = m + 1;
      end
    end
    m = m -1;
    nrsCyc = m;
    lrZC = prCyc(m);
  else
    nrsCyc = 1;
    lrZC = pZC(1);
    rsCyc = 0;
  end
  if (m ==1)
    nrsCyc = NaN;
    prCyc = NaN;
  end
end

function [drFirstZC, odrFirstZC, oodrFirstZC, drSecondZC, odrSecondZC, oodrSecondZC] = drZCFind(nrdFirstZC, nrdSecondZC, prFirstZC, prSecondZC)

  [drFirstZC(1), drSecondZC(1), odrFirstZC(1), odrSecondZC(1)] = AssignVal2Var(0);
  [oodrFirstZC(1), oodrSecondZC(1)] = AssignVal2Var(99999)
  [j,k,l] = AssignVal2Var(1);           % iterators

  while ((j+k) <= nrdFirstZC)&&(j<=nrdSecondZC)            % start at prFirstZC(2) (to check repetition)
    if prFirstZC(j+k)  <= prSecondZC(j)
      drFirstZC(k) = j + k;            % the ZC that is repeated
      odrFirstZC(k) = j + k - 1;           % the ZC before that
      oodrSecondZC(l) = j - 1;            % the ZC of the other Polarity before the repetition
      l = l + 1;
      k = k + 1;
    else
      j = j + 1;
    end
  end
  m = j + k;            % iterator from current position to end
  if m <= nrdFirstZC           % check for repeated nZC at the end.
    for m=m:nrdFirstZC
      drFirstZC(k) = m;
      odrFirstZC(k) = m - 1;
      oodrSecondZC(l) = 0;            % fill array with NULL because of later checks at those indices
      l = l + 1;
      k = k + 1;
    end
  end
  odrFirstZC(k) = 0;
  oodrSecondZC(l) = 0;
  oodrSecondZC(l+1) = 0;
  j = 1;
  k = 1;
  l = 1;

  while ((j+k) <= nrdSecondZC)&&(j<=nrdFirstZC)            % start at prSecondZC(2) (to check repetition)
    if prSecondZC(j+k)  <= prFirstZC(j)
      drSecondZC(k) = j + k;            % the ZC that is repeated
      odrSecondZC(k) = j + k - 1;           % the ZC before that
      oodrFirstZC(l) = j - 1;            % the ZC of the other Polarity before the repetition
      l = l + 1;
      k = k + 1;
    else
      j = j + 1;
    end
  end
  m = j + k;            % iterator from current position to end
  if m <= nrdSecondZC           % check for repeated nZC at the end.
    for m=m:nrdSecondZC
      drSecondZC(k) = m;
      odrSecondZC(k) = m - 1;
      oodrFirstZC(l) = 0;            % fill array with NULL because of later checks at those indices
      l = l + 1;
      k = k + 1;
    end
  end
  odrSecondZC(k) = 0;
  oodrFirstZC(l) = 0;
  oodrFirstZC(l+1) = 0;
end


function [N_PosHalfCyc, N_NegHalfCyc, N_Cyc, aN_PosHalfCyc, aN_NegHalfCyc, aN_Cyc, n_Cyc] = FindWholeCycs(frZC, nrPosCyc, nrNegCyc, prPosZC, prNegZC, odrPosZC, odrNegZC);
  aN_NegHalfCyc = 0;
  aN_PosHalfCyc = 0;
  aN_Cyc = 0;
  j = 1;
  k = 0;
  l = 0;
  m = 1;
  n = 1;

  if frZC == 0       % case: first ZC is negative --> Cycles get started at NegZCs
    while ((j+k+1) <= nrNegCyc)&&((j+l+1) <= nrPosCyc)
      if (j+k) == odrNegZC(k+1)
        k = k + 1;
      elseif (j+l) == odrPosZC(l+1)
        l = l + 1;
      else
        N_PosHalfCyc(j) = prPosZC(j+l) - prNegZC(j+k);
        N_NegHalfCyc(j) = prNegZC(j+k+1) - prPosZC(j+l);
        N_Cyc(j) = N_PosHalfCyc(j) + N_NegHalfCyc(j);
        aN_PosHalfCyc = aN_PosHalfCyc + N_PosHalfCyc(j);
        aN_NegHalfCyc = aN_NegHalfCyc + N_NegHalfCyc(j);
        aN_Cyc = aN_Cyc + N_Cyc(j);
        j = j + 1;
      end
    end
  else
    while ((j+k+1) <= nrPosCyc)&&((j+l+1) <= nrNegCyc)
      if (j+k) == odrPosZC(k+1)
        k = k + 1;
      elseif (j+l) == odrNegZC(l+1)
        l = l + 1;
      else
        N_NegHalfCyc(j) = prNegZC(j+l) - prPosZC(j+k);
        N_PosHalfCyc(j) = prPosZC(j+k+1) - prNegZC(j+l);
        N_Cyc(j) = N_NegHalfCyc(j) + N_PosHalfCyc(j);
        aN_PosHalfCyc = aN_PosHalfCyc + N_PosHalfCyc(j);
        aN_NegHalfCyc = aN_NegHalfCyc + N_NegHalfCyc(j);
        aN_Cyc = aN_Cyc + N_Cyc(j);
        j = j + 1;
      end
    end
  end
  n_Cyc = j - 1;
end

function [sCyc, asCyc, nCyc] = FindPolCyc(pZC, nZC);
  j = 1;
  asCyc = 0;
  for j=1:(nZC - 1)
    sCyc(j) = pZC(j + 1) - pZC(j);
    asCyc = asCyc + sCyc(j);
  end
  nCyc = j - 1;
end

function [nPeriod, Freq] = StatFilter(nPosCyc, nNegCyc)
  global fs;
  mnNegCyc = median(nNegCyc);
  mnPosCyc = median(nPosCyc);
  anNegCyc = mean(nNegCyc);
  anPosCyc = mean(nPosCyc);

  diffnNegCyc  = abs(mnNegCyc - anNegCyc);
  diffnPosCyc  = abs(mnPosCyc - anPosCyc);

  if (diffnPosCyc <= diffnNegCyc)
    nCycStable = (mnPosCyc + anPosCyc) / 2;
  else
    nCycStable = (mnNegCyc + anNegCyc) / 2;
  end

  OutlierArray(1) = abs(nCycStable - anPosCyc);
  OutlierArray(2) = abs(nCycStable - mnPosCyc);
  OutlierArray(3) = abs(nCycStable - anNegCyc);
  OutlierArray(4) = abs(nCycStable - mnNegCyc);
  [Maximum,Outlier] = max(OutlierArray);

  mArray(1) = anPosCyc;
  mArray(2) = mnPosCyc;
  mArray(3) = anNegCyc;
  mArray(4) = mnNegCyc;
  mArray([Outlier]) = [];           % remove outlier

  nPeriod = median(mArray);
  Freq = fs / nPeriod;
end

function [Freq, asCyc] = SinglePolStatfilter(sCyc)            % 'StatFilter' for only one polaroty 
  global fs;
  mnCyc = median(sCyc);
  anCyc = mean(sCyc);

  asCyc = (mnCyc + anCyc) / 2;
  Freq = fs / asCyc;
end

function [bReSizer, eReSizer, csWindow] = WinSize(cWindow, sWindow, sbSignal)
  bReSizer = ((cWindow - 1) * sWindow) + 1;
  eReSizer = cWindow * sWindow;
  if eReSizer > sbSignal           % case: last window (which is <csWindow)
    eReSizer = sbSignal;
  end
  csWindow = eReSizer - bReSizer + 1;
end

function [nMax, nMaxIndex] = kLargestElement(Array, n)
  nArray = Size(Array);
  TempArray = Array;
  for i=1:n
    [UpLimitIndex, UpLimit] = max(TempArray);
    if(i~=n)
      TempArray = 0;
    end
  end
  nMax = 0;
  nMaxIndex = 1;
  for i=1:nArray
    if (Array(i) < UpLimit) && (Array > nMax)
      nMax = Array(i);
      nMaxIndex = i;
    end
  end
end

function [PrintError] = PrintOutput(FileName, DataString, PrintError)
  global DataOutputDir;
  fid = fopen(strcat(DataOutputDir, FileName), 'w+t');
  if (fid ~= -1)
    fprintf(fid, DataString, '-ascii');
    fclose(fid);
    PrintError = PrintError;
  else
    PrintError = PrintError + 1;
  end
end


% Functions which shouldn't be necessary +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [SizeVar] = Size(Var)            % Matlabs size-function is automatically 2 dimensional.
  SizeVar = size(Var);
  SizeVar = SizeVar(2);
end

function [Var1, Var2, Var3, Var4, Var5, Var6, Var7, Var8, Var9] = InitString() % only single assignments per line in Matlab.
  Var1 = '';
  Var2 = '';
  Var3 = '';
  Var4 = '';
  Var5 = '';
  Var6 = '';
  Var7 = '';
  Var8 = '';
  Var9 = '';
end

function [Var1, Var2, Var3, Var4, Var5, Var6, Var7, Var8, Var9] = AssignVal2Var(Val) % only single assignments per line in Matlab.
  Var1 = Val;
  Var2 = Val;
  Var3 = Val;
  Var4 = Val;
  Var5 = Val;
  Var6 = Val;
  Var7 = Val;
  Var8 = Val;
  Var9 = Val;
end
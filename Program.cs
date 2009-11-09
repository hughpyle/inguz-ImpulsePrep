using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Globalization;
using DSPUtil;

// Copyright (c) 2006-2009 by Hugh Pyle, inguzaudio.com


// this code is particularly ugly, sorry to anyone trying to read it :-(
// there are lots of "quote"experimental"quote" and otherwise broken things

namespace ImpulsePrep
{
    /// <summary>
    /// ImpulsePrep
    /// commandline app for impulse-response preparation from sweep recordings
    /// and DRC for creation of the correction filters
    /// </summary>
    class Program
    {
        static char slash = Path.DirectorySeparatorChar;

        // config
        static int _nInFiles = 0;
        static string _thisFolder;
        static string _impulsesFolder;
        static List<string> _impulseFiles = new List<string>();
        static List<string> _filterFiles = new List<string>();
        static List<string> _tempFiles = new List<string>();
        static int _filterLen = 32768;
        static bool _noNorm = false;
        static bool _noDRC = false;
        static bool _noSkew = false;
        static bool _dbl = false;
        static bool _pcm = false;
        static bool _split = false;
        static bool _keepTempFiles = false;
        static bool _copy;
        static bool _env = false;
        static int _fmin = 0;
        static int _fmax = int.MaxValue;
        static bool _fminSpecified = false;
        static bool _fmaxSpecified = false;
        static bool _returnAll = false;
        static int _power = 0;
        static bool _doDirectFilters = false;
        static bool _noOverrideDRC = false;
        static bool _noSubsonicFilter = false;
        static int _peakPosL;
        static int _peakPosR;
        static int _refchannel = -1;
        static double _gain = Double.NaN; // dB.  NaN means normalize instead.

        // state
        static uint _sampleRate = 44100;


        static void Main(string[] args)
        {
            // Find where this executable is launched from
            string[] cargs = Environment.GetCommandLineArgs();
            _thisFolder = Path.GetDirectoryName(cargs[0]);
            if (String.IsNullOrEmpty(_thisFolder))
            {
                _thisFolder = Environment.CurrentDirectory;
            }

            string appData = Environment.GetFolderPath(Environment.SpecialFolder.CommonApplicationData);
            _impulsesFolder = Path.GetFullPath(Path.Combine(appData, "InguzEQ" + slash + "Impulses" + slash));

            string[] inFiles = new string[4];
            string inL = "";
            string inR = "";
            if (!DisplayInfo())
            {
                return;
            }

            bool ok = (args.Length > 0);
            bool longUsage = false;

            for (int j = 0; ok && j < args.Length; j++)
            {
                string arg = args[j];
                switch (args[j].ToUpperInvariant())
                {
                    case "/?":
                    case "-?":
                    case "/H":
                    case "/HELP":
                        ok = false;
                        longUsage = true;
                        break;

                    case "/L":
                    case "/0":
                        inFiles[0] = args[++j];
                        _nInFiles = Math.Max(_nInFiles, 1);
                        break;

                    case "/R":
                    case "/1":
                        inFiles[1] = args[++j];
                        _nInFiles = Math.Max(_nInFiles, 2);
                        break;

                    case "/2":
                        inFiles[2] = args[++j];
                        _nInFiles = Math.Max(_nInFiles, 3);
                        break;
                    case "/3":
                        inFiles[3] = args[++j];
                        _nInFiles = Math.Max(_nInFiles, 4);
                        break;

                    case "/LENGTH":
                        _filterLen = int.Parse(args[++j], CultureInfo.InvariantCulture);
                        if (_filterLen < 16)
                        {
                            throw new Exception("Length is too small.");
                        }
                        break;

                    case "/DBL":
                        _dbl = true;
                        break;

                    case "/PCM":
                        _pcm = true;
                        break;

                    case "/NODRC":
                        _noDRC = true;
                        break;

                    case "/NOSKEW":
                        _noSkew = true;
                        break;

                    case "/NONORM":
                        // No normalization of the impulse response (undocumented)
                        _noNorm = true;
                        break;

                    case "/SPLIT":
                        _split = true;
                        break;

                    case "/COPY":
                        _copy = true;
                        break;

                    case "/GAIN":
                        _gain = double.Parse(args[++j], CultureInfo.InvariantCulture);
                        break;

                    case "/ALL":
                        // Returns negative-time components as part of the impulse response
                        // (experimental, to be used for THD measurement)
                        _returnAll = true;
                        break;

                    case "/POWER":
                        // Raises sweep to power n
                        // (experimental, to be used for THD measurement)
                        _power = int.Parse(args[++j], CultureInfo.InvariantCulture);
                        break;

                    case "/FMIN":
                        // (experimental, i.e. broken)
                        _fmin = int.Parse(args[++j], CultureInfo.InvariantCulture);
                        _fminSpecified = true;
                        break;

                    case "/FMAX":
                        // (experimental, i.e. broken)
                        _fmax = int.Parse(args[++j], CultureInfo.InvariantCulture);
                        _fmaxSpecified = true;
                        break;

                    case "/DIRECT":
                        // Create filtered (direct-sound) filters
                        _doDirectFilters = true;
                        break;

                    case "/NOSUB":
                        // Don't apply subsonic filter to the impulse response
                        _noSubsonicFilter = true;
                        break;

                    case "/NOOVER":
                        // Don't override DRC's settings for filter type and length
                        _noOverrideDRC = true;
                        break;

                    case "/KEEPTEMP":
                        // Undocumented
                        _keepTempFiles = true;
                        break;

                    case "/REFCH":
                        // Override the reference-channel detection
                        _refchannel = int.Parse(args[++j], CultureInfo.InvariantCulture);
                        if (_refchannel<0 || _refchannel > _nInFiles - 1)
                        {
                            throw new Exception(String.Format("RefCh can only be from 0 to {0}.", _nInFiles-1));
                        }
                        break;

                    case "/ENV":
                        // Undocumented.  Save the Hilbert envelope
                        _env = true;
                        break;

                    case "-":
                        // ignore
                        break;

                    default:
                        ok = false;
                        break;
                }
            }
            if (!ok)
            {
                DisplayUsage(longUsage);
            }
            else
            {
                try
                {
                    if (!_noDRC)
                    {
                        if (!File.Exists(GetDRCExe()))
                        {
                            stderr.WriteLine("Denis Sbragion's DRC (http://drc-fir.sourceforge.net/) was not found.");
                            stderr.WriteLine("Only the impulse response will be calculated, not correction filters.");
                            stderr.WriteLine("");
                            _noDRC = true;
                        }
                    }
                    if (!_noDRC)
                    {
                        FileInfo[] drcfiles = new DirectoryInfo(_thisFolder).GetFiles("*.drc");
                        if (drcfiles.Length == 0)
                        {
                            stderr.WriteLine("No .drc files were found in the current folder.");
                            stderr.WriteLine("Only the impulse response will be calculated, not correction filters.");
                            stderr.WriteLine("");
                            _noDRC = true;
                        }
                    }


                    for(int i=0; i<_nInFiles; i++)
                    {
                        string inFile = inFiles[i];
                        if (String.IsNullOrEmpty(inFile))
                        {
                            stderr.WriteLine("Error: The {0} input file was not specified.", FileDescription(i));
                            return;
                        }
                        if (!File.Exists(inFile))
                        {
                            stderr.WriteLine("Error: The {0} input file {1} was not found.", FileDescription(i), inFile);
                            return;
                        }

                        for (int j = 0; j < i; j++)
                        {
                            if (inFile.Equals(inFiles[j]))
                            {
                                stderr.WriteLine("Warning: The same input file ({0}) was specified for both {1} and {2}!", inFile, FileDescription(j), FileDescription(i));
                                //stderr.WriteLine();
                            }
                        }
                    }


                    // Temporary
                    if (_nInFiles != 2)
                    {
                        stderr.WriteLine("Error: Two input files must be specified.");
                        return;
                    }
                    inL = inFiles[0];
                    inR = inFiles[1];
                    // end temporary


                    uint sampleRate;
                    List<SoundObj> impulses;
                    List<ISoundObj> filteredImpulses;
                    List<string> impDirects;
                    List<Complex[]> impulseFFTs;
                    List<double> maxs;

                    SoundObj impulseL;
                    SoundObj impulseR;
                    ISoundObj filteredImpulseL = null;
                    ISoundObj filteredImpulseR = null;

                    string impDirectL = null;
                    string impDirectR = null;

                    Complex[] impulseLFFT;
                    Complex[] impulseRFFT;
                    WaveWriter writer;

                    ISoundObj buff;
                    double g;

                    if (!_keepTempFiles)
                    {
                        _tempFiles.Add("rps.pcm");
                        _tempFiles.Add("rtc.pcm");
                    }

                    // Find the left impulse
                    stderr.WriteLine("Processing left measurement ({0})...", inL);
                    impulseL = Deconvolve(inL, out impulseLFFT, out _peakPosL);
                    sampleRate = impulseL.SampleRate;
                    _sampleRate = sampleRate;
                    double peakM = Math.Round(MathUtil.Metres(_peakPosL, sampleRate), 2);
                    double peakFt = Math.Round(MathUtil.Feet(_peakPosL, sampleRate), 2);
                    stderr.WriteLine("  Impulse peak at sample {0} ({1}m, {2}ft)", _peakPosL, peakM, peakFt);

                    // Write to PCM
                    string impFileL = Path.GetFileNameWithoutExtension(inL) + "_imp" + ".pcm";
                    if (!_keepTempFiles)
                    {
                        _tempFiles.Add(impFileL);
                    }
                    writer = new WaveWriter(impFileL);
                    writer.Input = impulseL;
                    writer.Format = WaveFormat.IEEE_FLOAT;
                    writer.BitsPerSample = 32;
                    writer.SampleRate = _sampleRate;
                    writer.Raw = true;
                    writer.Run();
                    writer.Close();


                    // Write the impulseFFT to disk
                    int L = impulseLFFT.Length;
                    string impTempL = Path.GetFileNameWithoutExtension(inL) + "_imp" + ".dat";
                    _tempFiles.Add(impTempL);
                    writer = new WaveWriter(impTempL);
                    writer.Input = new CallbackSource(2, sampleRate, delegate(long j)
                    {
                        if (j >= L / 2)
                        {
                            return null;
                        }
                        Complex si = impulseLFFT[j]; // +impulseLFFT[L - j - 1];
                        ISample s = new Sample2();
                        s[0] = si.Magnitude;
                        s[1] = si.Phase / Math.PI;
                        return s;
                    });
                    writer.Format = WaveFormat.IEEE_FLOAT;
                    writer.BitsPerSample = 32;
                    writer.SampleRate = _sampleRate;
                    writer.Raw = false;
                    writer.Run();
                    writer.Close();
                    writer = null;

                    impulseLFFT = null;
                    GC.Collect();

                    if (_doDirectFilters)
                    {
                        // Sliding low-pass filter over the impulse
                        stderr.WriteLine("  Filtering...");
                        filteredImpulseL = SlidingLowPass(impulseL, _peakPosL);

                        // Write PCM for the filtered impulse
                        impDirectL = Path.GetFileNameWithoutExtension(inL) + "_impfilt" + ".pcm";
                        if (!_keepTempFiles)
                        {
                            _tempFiles.Add(impDirectL);
                        }
                        writer = new WaveWriter(impDirectL);
                        writer.Input = filteredImpulseL;
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.SampleRate = _sampleRate;
                        writer.BitsPerSample = 32;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        writer = null;
                        filteredImpulseL.Reset();
                    }

                    GC.Collect();
                    stderr.WriteLine("  Deconvolution: left impulse done.");
                    stderr.WriteLine();

                    // Find the right impulse
                    stderr.WriteLine("Processing right measurement ({0})...", inR);
                    impulseR = Deconvolve(inR, out impulseRFFT, out _peakPosR);
                    peakM = Math.Round(MathUtil.Metres(_peakPosR, sampleRate), 2);
                    peakFt = Math.Round(MathUtil.Feet(_peakPosR, sampleRate), 2);
                    stderr.WriteLine("  Impulse peak at sample {0} ({1}m, {2}ft)", _peakPosR, peakM, peakFt);

                    // Write to PCM
                    string impFileR = Path.GetFileNameWithoutExtension(inR) + "_imp" + ".pcm";
                    if (!_keepTempFiles)
                    {
                        _tempFiles.Add(impFileR);
                    }
                    writer = new WaveWriter(impFileR);
                    writer.Input = impulseR;
                    writer.Format = WaveFormat.IEEE_FLOAT;
                    writer.BitsPerSample = 32;
                    writer.SampleRate = _sampleRate;
                    writer.Raw = true;
                    writer.Run();
                    writer.Close();

                    
                    // Write the impulseFFT magnitude to disk
                    L = impulseRFFT.Length;
                    string impTempR = Path.GetFileNameWithoutExtension(inR) + "_imp" + ".dat";
                    _tempFiles.Add(impTempR);
                    writer = new WaveWriter(impTempR);
                    writer.Input = new CallbackSource(2, impulseR.SampleRate, delegate(long j)
                    {
                        if (j >= L / 2)
                        {
                            return null;
                        }
                        Complex si = impulseRFFT[j]; // +impulseRFFT[L - j - 1];
                        ISample s = new Sample2();
                        s[0] = si.Magnitude;
                        s[1] = si.Phase / Math.PI;
                        return s;
                    });
                    writer.Format = WaveFormat.IEEE_FLOAT;
                    writer.BitsPerSample = 32;
                    writer.SampleRate = _sampleRate;
                    writer.Raw = false;
                    writer.Run();
                    writer.Close();
                    writer = null;

                    impulseRFFT = null;
                    GC.Collect();

                    if (_doDirectFilters)
                    {
                        // Sliding low-pass filter over the impulse
                        stderr.WriteLine("  Filtering...");
                        filteredImpulseR = SlidingLowPass(impulseR, _peakPosR);

                        // Write PCM for the filtered impulse
                        impDirectR = Path.GetFileNameWithoutExtension(inR) + "_impfilt" + ".pcm";
                        if (!_keepTempFiles)
                        {
                            _tempFiles.Add(impDirectR);
                        }
                        writer = new WaveWriter(impDirectR);
                        writer.Input = filteredImpulseR;
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        writer = null;
                        filteredImpulseR.Reset();
                    }

                    GC.Collect();

                    stderr.WriteLine("  Deconvolution: right impulse done.");
                    stderr.WriteLine();

                    // Join the left and right impulse files (truncated at 65536) into a WAV
                    // and normalize loudness for each channel
                    stderr.WriteLine("Splicing and normalizing (1)");
                    ChannelSplicer longstereoImpulse = new ChannelSplicer();

                    // (Don't normalize each channel's volume separately if _returnAll, it's just too expensive)
                    if (_returnAll)
                    {
                        buff = impulseL;
                    }
                    else
                    {
                        buff = new SoundBuffer(new SampleBuffer(impulseL).Subset(0, 131071));
                        g = Loudness.WeightedVolume(buff);
                        (buff as SoundBuffer).ApplyGain(1 / g);
                    }
                    longstereoImpulse.Add(buff);

                    if (_returnAll)
                    {
                        buff = impulseR;
                    }
                    else
                    {
                        buff = new SoundBuffer(new SampleBuffer(impulseR).Subset(0, 131071));
                        g = Loudness.WeightedVolume(buff);
                        (buff as SoundBuffer).ApplyGain(1 / g);
                    }
                    longstereoImpulse.Add(buff);

                    ISoundObj stereoImpulse = longstereoImpulse;

                    _impulseFiles.Add("Impulse_Response_Measured.wav: stereo impulse response from measurements");
                    writer = new WaveWriter("Impulse_Response_Measured.wav");
                    writer.Input = longstereoImpulse;
                    writer.Format = WaveFormat.IEEE_FLOAT;
                    writer.BitsPerSample = 32;
                    writer.SampleRate = _sampleRate;
                    writer.Normalization = -1;
                    writer.Raw = false;
                    writer.Run();
                    writer.Close();
                    writer = null;

                    if (_env)
                    {
                        // Also save the Hilbert envelope
                        HilbertEnvelope env = new HilbertEnvelope(8191);
                        env.Input = longstereoImpulse;
                        _impulseFiles.Add("Impulse_Response_Envelope.wav: Hilbert envelope of the impulse response");
                        writer = new WaveWriter("Impulse_Response_Envelope.wav");
                        writer.Input = env;
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Normalization = -1;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        writer = null;
                    }
                    
                    if (_dbl)
                    {
                        // Create DBL files for Acourate
                        _impulseFiles.Add("PulseL.dbl: impulse response, raw data (64-bit float), left channel ");
                        _impulseFiles.Add("PulseR.dbl: impulse response, raw data (64-bit float), right channel");
                        _impulseFiles.Add("  (use skew=" + (_peakPosL - _peakPosR) + " for time alignment)");
                        WriteImpulseDBL(stereoImpulse, "PulseL.dbl", "PulseR.dbl");
                    }

                    if (_pcm)
                    {
                        // Create PCM files for Octave (etc)
                        _impulseFiles.Add("LUncorrected.pcm: impulse response, raw data (32-bit float), left channel");
                        _impulseFiles.Add("RUncorrected.pcm: impulse response, raw data (32-bit float), right channel");
                        WriteImpulsePCM(stereoImpulse, "LUncorrected.pcm", "RUncorrected.pcm");
                    }

                    stereoImpulse = null;
                    longstereoImpulse = null;
                    buff = null;
                    GC.Collect();

                    if (_doDirectFilters)
                    {
                        // Same for the filtered impulse response
                        stderr.WriteLine("Splicing and normalizing (2)");
                        ChannelSplicer longstereoImpulseF = new ChannelSplicer();

                        buff = new SoundBuffer(new SampleBuffer(filteredImpulseL).Subset(0, 131071));
                        double gL = Loudness.WeightedVolume(buff);
                        (buff as SoundBuffer).ApplyGain(1 / gL);
                        longstereoImpulseF.Add(buff);
                        FilterProfile lfgDirectL = new FilterProfile(buff, 0.5);

                        buff = new SoundBuffer(new SampleBuffer(filteredImpulseR).Subset(0, 131071));
                        double gR = Loudness.WeightedVolume(buff);
                        (buff as SoundBuffer).ApplyGain(1 / gR);
                        longstereoImpulseF.Add(buff);
                        FilterProfile lfgDirectR = new FilterProfile(buff, 0.5);

                        _impulseFiles.Add("Impulse_Response_Filtered.wav: approximation to direct-sound impulse response");
                        writer = new WaveWriter("Impulse_Response_Filtered.wav");
                        writer.Input = longstereoImpulseF;
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Normalization = -1;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        double gg = writer.Gain;
                        writer = null;
                        longstereoImpulseF = null;


                        ChannelSplicer longstereoImpulseD = new ChannelSplicer();

                        Mixer diffuse = new Mixer();
                        diffuse.Add(impulseL, 1.0);
                        diffuse.Add(filteredImpulseL, -1.0);
                        buff = new SoundBuffer(new SampleBuffer(diffuse).Subset(0, 131071));
                        (buff as SoundBuffer).ApplyGain(1 / gL);
                        longstereoImpulseD.Add(buff);
                        FilterProfile lfgDiffuseL = new FilterProfile(buff, 0.5);

                        diffuse = new Mixer();
                        diffuse.Add(impulseR, 1.0);
                        diffuse.Add(filteredImpulseR, -1.0);
                        buff = new SoundBuffer(new SampleBuffer(diffuse).Subset(0, 131071));
                        (buff as SoundBuffer).ApplyGain(1 / gR);
                        longstereoImpulseD.Add(buff);
                        FilterProfile lfgDiffuseR = new FilterProfile(buff, 0.5);

                        _impulseFiles.Add("Impulse_Response_Diffuse.wav: approximation to diffuse-field remnant");
                        writer = new WaveWriter("Impulse_Response_Diffuse.wav");
                        writer.Input = longstereoImpulseD;
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Gain = gg;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        writer = null;

                        // Filter the diffuse-field curve against double the diffuse-field curve
                        FilterImpulse fiDiffuse = new FilterImpulse(8192, HRTF.diffuseDiff0() * 2, FilterInterpolation.COSINE, sampleRate);
                        FastConvolver co = new FastConvolver(longstereoImpulseD, fiDiffuse);
                        SoundBuffer buffd = new SoundBuffer(co);
                        _impulseFiles.Add("Impulse_Response_Diffuse_Comp.wav: filtered diffuse-field remnant");
                        writer = new WaveWriter("Impulse_Response_Diffuse_Comp.wav");
                        writer.Input = buffd.Subset(4096);
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Gain = gg;
                        writer.Raw = false;
                        writer.Run();
                        writer.Close();
                        writer = null;

                        longstereoImpulseD = null;

                        bool any = false;
                        string jsonFile = "Diff.json";
                        FileStream fs = new FileStream(jsonFile, FileMode.Create);
                        StreamWriter sw = new StreamWriter(fs);
                        sw.WriteLine("{");
                        FilterProfile lfgDiffL = lfgDirectL - lfgDiffuseL;
                        if (lfgDiffL != null)
                        {
                            if (any) sw.WriteLine(",");
                            any = true;
                            sw.Write(lfgDiffL.ToJSONString("DiffL", "Diffuse field relative to direct, left channel"));
                        }
                        FilterProfile lfgDiffR = lfgDirectR - lfgDiffuseR;
                        if (lfgDiffR != null)
                        {
                            if (any) sw.WriteLine(",");
                            any = true;
                            sw.Write(lfgDiffR.ToJSONString("DiffR", "Diffuse field relative to direct, right channel"));
                        }
                        sw.WriteLine("}");
                        sw.Close();
                        fs.Close();
                    }
                    buff = null;
                    GC.Collect();

                    System.Console.Error.WriteLine();

                    if (!_noDRC)
                    {
                        // Analyze the freq response
                        // and create targets
                        // target_full.txt and target_half.txt
                        stderr.WriteLine("Analyzing response curves.");
                        Prep(impTempL, impTempR, "Impulse_Response_Measured.wav", "NoCorrection");

                        // Call DRC to create the filters
                        // then splice the DRC left & right output files together
                        stderr.WriteLine("Preparing for DRC.");
                        if (DoDRC(impFileL, impFileR, impDirectL, impDirectR, _peakPosL, _peakPosR, "Impulse_Response_Measured.wav", "Impulse_Response_Filtered.wav"))
                        {
                            stderr.WriteLine("Success!");
                        }
                    }

                    // Report names of the impulse files created
                    if (_impulseFiles.Count == 0)
                    {
                        System.Console.Error.WriteLine("No impulse response files were created.");
                    }
                    if (_impulseFiles.Count > 0)
                    {
                        System.Console.Error.WriteLine("Impulse response files were created:");
                        foreach (string f in _impulseFiles)
                        {
                            string s = "  " + f;
                            System.Console.Error.WriteLine(s);
                        }
                    }


                    // Report names of the filter files created
                    if (_filterFiles.Count == 0 && !_noDRC)
                    {
                        System.Console.Error.WriteLine("No correction filter files were created.");
                    }
                    if (_filterFiles.Count > 0)
                    {
                        System.Console.Error.WriteLine("Correction filter files were created:");
                        foreach (string f in _filterFiles)
                        {
                            string s = "  " + f;
                            if (_copy)
                            {
                                try
                                {
                                    File.Copy(f, Path.Combine(_impulsesFolder, f), true);
                                    s += " (copied)";
                                }
                                catch (Exception e)
                                {
                                    s += " (not copied: " + e.Message + ")";
                                }
                            }
                            System.Console.Error.WriteLine(s);
                        }
                    }
                    if (_peakPosL == _peakPosR)
                    {
                        System.Console.Error.WriteLine();
                        System.Console.Error.WriteLine("Zero time difference between channels.  Are you sure the recordings are correct?");
                    }
                }
                catch (Exception e)
                {
                    stderr.WriteLine();
                    stderr.WriteLine(e.Message);
                    stderr.WriteLine(e.StackTrace);
                }
                finally
                {
                    foreach (string tempFile in _tempFiles)
                    {
                        try
                        {
                            File.Delete(tempFile);
                        }
                        catch (Exception) { /* ignore */ }
                    }
                }
            }
            stderr.Flush();
        }

        static string FileDescription(int nFile)
        {
            if (_nInFiles == 2)
            {
                return nFile == 0 ? "left" : "right";
            }
            else
            {
                return String.Format("channel {0}", nFile);
            }
        }

        static void WriteImpulseDBL(ISoundObj stereoImpulse, string fileNameL, string fileNameR)
        {
            ISoundObj sb;
            WaveWriter writer;
            ushort bits = 64;
            int nTot = 65535;
            int nBef = 6000;

            // need padding before the sample.
            // window the sample first, then pad it.
            // Total length 65536, includes nPad and _peakPosL before impulse
            int nPad = nBef - _peakPosL;
            int usableLength = nTot - nPad;
            int windowCenter = usableLength / 2;
            int windowSide = _peakPosL * 2 / 3;
            int windowFlat = windowCenter - windowSide;
            ISoundObj ch = new SingleChannel(stereoImpulse, 0);
            ISoundObj wn = new BlackmanHarris(windowCenter, windowSide, windowFlat);
            wn.Input = ch;
            sb = new SampleBuffer(wn).PaddedSubset(-nPad, nTot);

            writer = new WaveWriter(fileNameL);
            writer.Input = sb;
            writer.Format = WaveFormat.IEEE_FLOAT;
            writer.BitsPerSample = bits;
            writer.SampleRate = _sampleRate;
            writer.Raw = true;
            writer.Run();
            writer.Close();
            writer = null;

            nPad = nBef - _peakPosR;
            usableLength = nTot - nPad;
            windowCenter = usableLength / 2;
            windowSide = _peakPosR * 2 / 3;
            windowFlat = windowCenter - windowSide;
            ch = new SingleChannel(stereoImpulse, 1);
            wn = new BlackmanHarris(windowCenter, windowSide, windowFlat);
            wn.Input = ch;
            sb = new SampleBuffer(wn).PaddedSubset(-nPad, nTot);

            writer = new WaveWriter(fileNameR);
            writer.Input = sb;
            writer.Format = WaveFormat.IEEE_FLOAT;
            writer.BitsPerSample = bits;
            writer.SampleRate = _sampleRate;
            writer.Raw = true;
            writer.Run();
            writer.Close();
            writer = null;
        }

        static void WriteImpulsePCM(ISoundObj stereoImpulse, string fileNameL, string fileNameR)
        {
            ISoundObj sb;
            WaveWriter writer;
            ushort bits = 32;
            int nTot = 196607;
            int nBef = 98300;

            // need padding before the sample.
            // window the sample first, then pad it.
            // Total length 65536, includes nPad and _peakPosL before impulse
            int nPad = nBef - _peakPosL;
            int usableLength = nTot - nPad;
            int windowCenter = usableLength / 2;
            int windowSide = _peakPosL * 2 / 3;
            int windowFlat = windowCenter - windowSide;
            ISoundObj ch = new SingleChannel(stereoImpulse, 0);
            ISoundObj wn = new BlackmanHarris(windowCenter, windowSide, windowFlat);
            wn.Input = ch;
            sb = new SampleBuffer(wn).PaddedSubset(-nPad, nTot);

            writer = new WaveWriter(fileNameL);
            writer.Input = sb;
            writer.Format = WaveFormat.IEEE_FLOAT;
            writer.BitsPerSample = bits;
            writer.SampleRate = _sampleRate;
            writer.Raw = true;

            writer.Normalization = -6.0;
            writer.NormalizationType = NormalizationType.PEAK_DBFS;

            writer.Run();
            writer.Close();

            double g = writer.Gain;

            writer = null;

            nPad = nBef - _peakPosR;
            usableLength = nTot - nPad;
            windowCenter = usableLength / 2;
            windowSide = _peakPosR * 2 / 3;
            windowFlat = windowCenter - windowSide;
            ch = new SingleChannel(stereoImpulse, 1);
            wn = new BlackmanHarris(windowCenter, windowSide, windowFlat);
            wn.Input = ch;
            sb = new SampleBuffer(wn).PaddedSubset(-nPad, nTot);

            writer = new WaveWriter(fileNameR);
            writer.Input = sb;
            writer.Format = WaveFormat.IEEE_FLOAT;
            writer.BitsPerSample = bits;
            writer.SampleRate = _sampleRate;
            writer.Raw = true;
            writer.Gain = g;
            writer.Run();
            writer.Close();
            writer = null;
        }

        static ISoundObj LowPassFiltered(ISoundObj input, double freqStart, double gainEnd)
        {
            uint sampleRate = input.SampleRate;
            FilterProfile lfg = new FilterProfile();
            lfg.Add(new FreqGain(freqStart, 0));
            lfg.Add(new FreqGain(0.499*sampleRate, gainEnd));
            FastConvolver conv = new FastConvolver();
            conv.Input = input;
            conv.impulse = new FilterImpulse(8192, lfg, FilterInterpolation.COSINE, sampleRate);
            return conv;
        }

        static ISoundObj SlidingLowPass(ISoundObj impulse, int peakPos)
        {
            // Filter the impulse response with a sliding lowpass filter
            // - Take the first (peakpos) samples unchanged
            // - Take the next 2.5ms (~110 samples @ 44100) unchanged
            // - Fade to a lowpass-filtered version

            uint sampleRate = impulse.SampleRate;
            int startpos1 = (int)(sampleRate * 0.0025);     // 2.5mS, 110 samples
            int startpos2 = (int)(sampleRate * 0.0050);
            int startpos3 = (int)(sampleRate * 0.0100);
            int startpos4 = (int)(sampleRate * 0.0200);
            int startpos5 = (int)(sampleRate * 0.0400);

            List<IEnumerator<ISample>> filters = null;
            int nFilter = 0;
            int nFilters = 0;
            int x = 0;
            int f0len = 0;
            ISoundObj filteredImpulse = new CallbackSource(1, sampleRate, delegate(long j)
            {
                if(j==0)
                {
                    // initialize the list of filters
                    impulse.Reset();
                    filters = new List<IEnumerator<ISample>>(6);
                    ISoundObj f0 = LowPassFiltered(impulse, 20, 0);
                    f0len = f0.Iterations;
                    filters.Add(f0.Samples);        // unfiltered
                    filters.Add(LowPassFiltered(impulse, 640, -10).Samples);      // LPF from 640Hz, kick in at +110
                    filters.Add(LowPassFiltered(impulse, 320, -20).Samples);       // LPF from 320Hz, kick in at +220 (etc)
                    filters.Add(LowPassFiltered(impulse, 160, -30).Samples);       // +440
                    filters.Add(LowPassFiltered(impulse, 80, -40).Samples);        // +880
                    filters.Add(LowPassFiltered(impulse, 40, -50).Samples);        // +1760 = 40ms
                    nFilter = 0;
                    nFilters = filters.Count;
                }

                // Move the filters along a notch.
                // (or, right at the beginning, move all the filters along by a lot, compensating for their center)
                // For perf, we don't need to keep enumerating the filters we've finished using.
                bool moreData = true;
                int nSkip = 1;
                if (j == 0)
                {
                    nSkip = (f0len/2)+1;
                }
                for (int skip = 0; skip < nSkip; skip++)
                {
                    int nStart = (nFilter == 0) ? 0 : nFilter - 1;    // Keep moving even the old filters along so mixing works later
                    for (int f = nStart; f < nFilters; f++)
                    {
                        moreData &= filters[f].MoveNext();
                    }
                }
                if (!moreData) return null;

                // Decide which filter we'll use
                x++;
                if (j == peakPos + startpos1) { nFilter = 1; x = 0; }
                if (j == peakPos + startpos2) { nFilter = 2; x = 0; }
                if (j == peakPos + startpos3) { nFilter = 3; x = 0; }
                if (j == peakPos + startpos4) { nFilter = 4; x = 0; }
                if (j == peakPos + startpos5) { nFilter = 5; x = 0; }

                // Pick the sample from the current filter
                ISample s = filters[nFilter].Current;

                // For a short region after switch-over, mix slowly with the sample from the previous filter
                int w = 20;
                if (x < w && nFilter > 0)
                {
                    ISample sPrev = filters[nFilter-1].Current;
                    for (int c = 0; c < s.NumChannels; c++)
                    {
                        s[c] = ((double)x/w)*s[c] + ((double)(w-x)/w)*sPrev[c];
                    }
                }
                return s;
            });
            return filteredImpulse;
        }

        static FilterProfile bandpass(uint centerFreq, uint sr)
        {
            FilterProfile lfg = new FilterProfile();
            lfg.Add(new FreqGain(centerFreq * 0.5, -190));
            lfg.Add(new FreqGain(centerFreq * 1.00, 0));
            lfg.Add(new FreqGain(Math.Min(centerFreq * 2,sr), -190));
            return lfg;
        }

        static void FindPeaks(ISoundObj impulse)
        {
            // Input: a single-channel impulse

            // Find the peak positions:
            // - unfiltered
            // - filtered with various bandpass filters

            uint sr = impulse.SampleRate;
            ushort nc = impulse.NumChannels;

            double peakM = Math.Round(MathUtil.Metres(_peakPosL, sr), 2);
            double peakFt = Math.Round(MathUtil.Feet(_peakPosL, sr), 2);
            stderr.WriteLine("  Impulse peak at sample {0} ({1}m, {2}ft)", _peakPosL, peakM, peakFt);

            FilterImpulse fi;
            WaveWriter wri;
            FastConvolver co;

            fi = new FilterImpulse(2048, bandpass(400,sr), FilterInterpolation.COSINE, sr);
            co = new FastConvolver(impulse, fi);
            wri = new WaveWriter("bp_400.wav", nc, sr, 16, DitherType.NONE, WaveFormat.PCM);
            wri.Input = co;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();

            fi = new FilterImpulse(2048, bandpass(6000, sr), FilterInterpolation.COSINE, sr);
            co = new FastConvolver(impulse, fi);
            wri = new WaveWriter("bp_6000.wav", nc, sr, 16, DitherType.NONE, WaveFormat.PCM);
            wri.Input = co;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();

            // and

            fi = new FilterImpulse(2048, bandpass(160, sr), FilterInterpolation.COSINE, sr);
            co = new FastConvolver(impulse, fi);
            wri = new WaveWriter("bp_160.wav", nc, sr, 16, DitherType.NONE, WaveFormat.PCM);
            wri.Input = co;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();

            fi = new FilterImpulse(2048, bandpass(2560, sr), FilterInterpolation.COSINE, sr);
            co = new FastConvolver(impulse, fi);
            wri = new WaveWriter("bp_2560.wav", nc, sr, 16, DitherType.NONE, WaveFormat.PCM);
            wri.Input = co;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();

            fi = new FilterImpulse(2048, bandpass(18000, sr), FilterInterpolation.COSINE, sr);
            co = new FastConvolver(impulse, fi);
            wri = new WaveWriter("bp_18k.wav", nc, sr, 16, DitherType.NONE, WaveFormat.PCM);
            wri.Input = co;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();

        }

        static SoundObj Deconvolve(string infile, out Complex[] impulseFFT, out int peakpos)
        {
            WaveReader reader = new WaveReader(infile);
            ushort nChannels = reader.NumChannels;
            uint sampleRate = reader.SampleRate;

            CallbackSource cs;
            SingleChannel[] channels = new SingleChannel[2];
            Complex[][][] data = new Complex[2][][];
            double[] stdDev = new double[2];
            double[] maxLogs = new double[2];
            double[] maxs = new double[2];
            double[] avgLog = new double[2];

            Complex[] swp_data = null;
            Complex[] mea_data = null;

            if (_fmax == int.MaxValue)
            {
                _fmax = (int)sampleRate / 2;
            }
            bool isEstimatedSweepRange = false;

            if (nChannels != 2)
            {
                throw new Exception("Input must have two channels.");
            }

            peakpos = 0;
            int cSwp = 0;
            int cMea = 1;
            int L = 0;
            int Nh = 0;
            double mx;

            double max = 0;
            double maxLog = 0;
            double stdev = 0;
            double avg = 0;

            for (int iteration = 1; iteration <= 2; iteration++)
            {
                // Read and FFT all the data
                // one channel at a time
                for (ushort c = 0; c < 2; c++)
                {
                    SingleChannel channel = reader.Channel(c);
                    Complex[] cdata;

                    SoundBuffer buff = new SoundBuffer(channel);
                    buff.ReadAll();

                    if (iteration==2 && _split)
                    {
                        // Split up the input file
                        string infile2 = Path.ChangeExtension(Path.GetFileName(infile), ".PCM");
                        if (c == cSwp)
                        {
                            WaveWriter wri = new WaveWriter("refchannel_" + infile2, 1, channel.SampleRate, 32, DitherType.NONE, WaveFormat.IEEE_FLOAT);
                            wri.Input = buff;
                            wri.Run(); wri.Close();
                        }
                        if (c == cMea)
                        {
                            WaveWriter wri = new WaveWriter("sweep_" + infile2, 1, channel.SampleRate, 32, DitherType.NONE, WaveFormat.IEEE_FLOAT);
                            wri.Input = buff;
                            wri.Run(); wri.Close();
                        }
                    }

                    // And then double in length to prevent wraparound
                    buff.PadTo(buff.Count * 2);

                    // Pad to next higher power of two
                    buff.PadToPowerOfTwo();

                    // Read out into array of complex
                    data[c] = buff.ToComplexArray();

                    // Then we're done with the buffer for this channel
                    buff = null;
                    GC.Collect();

                    cdata = data[c][0];

                    if (iteration==2 && c==cSwp && _power > 0)
                    {
                        // Deconvolve against a power of the sweep,
                        // for distortion measurement of harmonic _power
                        Complex p = new Complex((double)_power,0);
                        for (int j = 0; j < cdata.Length; j++)
                        {
                            cdata[j].Pow(p);
                        }
                    }

                    // FFT in place
                    Fourier.FFT(cdata.Length, cdata);


                    if (false && iteration==1)
                    {
                        // write the fft magnitudes to disk
                        cs = new CallbackSource(1, sampleRate, delegate(long j)
                        {
                            if (j >= cdata.Length)
                            {
                                return null;
                            }
                            Complex si = cdata[j];
                            Sample s = new Sample(1);
                            double f = (double)j * sampleRate / cdata.Length;
                            s[0] = mag(sampleRate, f, si.Magnitude);
                            return s;
                        });
                        // cs.SampleRate = sampleRate;
                        // cs.NumChannels = 1;
                        WaveWriter writer = new WaveWriter("fft_" + c + "_" + infile);
                        writer.Format = WaveFormat.IEEE_FLOAT;
                        writer.BitsPerSample = 32;
                        writer.SampleRate = _sampleRate;
                        writer.Input = cs;
                        writer.Normalization = -3;
                        writer.Run();
                        writer.Close();
                    }


                    // Take a slice of the FFT, std dev of log(|fft|),
                    // the lower value should be the sweep
                    int n3 = cdata.Length / 4;
                    int n1 = n3 / 2;
                    int n2 = n1 + n3;
                    get_stddev(sampleRate, n1, n2, cdata, out max, out maxLog, out stdev, out avg);

                    maxs[c] = max;
                    maxLogs[c] = maxLog;
                    stdDev[c] = stdev;
                    avgLog[c] = avg;

                    // Trace.WriteLine("Channel {0} standard deviation {1}", c, stdDev[c]);
                }
                GC.Collect();

                if (iteration == 1)
                {
                    if (stdDev[1] < stdDev[0])
                    {
                        if (_refchannel == -1)
                        {
                            cSwp = 1;
                            cMea = 0;
                            stderr.WriteLine("  Right channel seems to be the sweep");
                        }
                        else
                        {
                            stderr.WriteLine("  Right channel seems to be the sweep");
                            stderr.WriteLine("  But you said use refchannel {0}, so using that.", _refchannel);
                            cSwp = _refchannel;
                            cMea = (_nInFiles - 1) - _refchannel;
                        }
                    }
                    else
                    {
                        if (_refchannel == -1)
                        {
                            stderr.WriteLine("  Left channel seems to be the sweep");
                        }
                        else
                        {
                            stderr.WriteLine("  Left channel seems to be the sweep");
                            stderr.WriteLine("  But you said use refchannel {0}, so using that.", _refchannel);
                            cSwp = _refchannel;
                            cMea = (_nInFiles - 1) - _refchannel;
                        }
                    }
                }

                swp_data = data[cSwp][0];
                mea_data = data[cMea][0];

                L = swp_data.Length;
                Nh = L / 2;

                // stderr.WriteLine("avgLog=" + avgLog[cSwp] + " maxLog=" + maxLogs[cSwp]);
                
                double hz1 = L / sampleRate;
                if (false && _fmin == 0)
                {
                    isEstimatedSweepRange = true;
                    // Working back from 100Hz,
                    // Look for the first range where average of a 1Hz window
                    // is less than 0.7* average for the sweep itself
                    int kmin = (int)(hz1 * 100);
                    _fmin = 100;
                    while (kmin > 0)
                    {
                        get_stddev(sampleRate, kmin, (int)(kmin + hz1), swp_data, out max, out maxLog, out stdev, out avg);
                        if (avg < avgLog[cSwp] * 0.5)
                        {
                            break;
                        }
                        kmin -= (int)hz1;
                        _fmin--;
                    }
                    // stderr.WriteLine("avg/2: kmin=" + kmin + ", _fmin=" + _fmin);
                }

                if (false && _fmax == sampleRate / 2)
                {
                    isEstimatedSweepRange = true;
                    // Working forward from (say) 15kHz,
                    // Look for the first range where average of a 100Hz window
                    // is less than 0.7* average for the sweep itself
                    int kmax = (int)(hz1 * 10000);
                    _fmax = 10000;
                    get_stddev(sampleRate, kmax, (int)(kmax + 100 * hz1), swp_data, out max, out maxLog, out stdev, out avg);
                    double stdTest = stdev;
                    while (kmax < L / 2)
                    {
                        get_stddev(sampleRate, kmax, (int)(kmax + 100 * hz1), swp_data, out max, out maxLog, out stdev, out avg);
                        if (avg < avgLog[cSwp] * 0.5)
                        {
                            break;
                        }
                        kmax += (int)(100 * hz1);
                        _fmax += 100;
                    }
                    // stderr.WriteLine("StdDev Peak: kmax=" + kmax + ", _fmax=" + _fmax + ", sdev=" + stdev + " ref " + stdTest + " (1Hz=" + hz1 + "), avgLog=" + avg);
                }

                if (!_noSubsonicFilter)
                {
                    // Filter LF from the measurement data
                    // to avoid spurious stuff below 15Hz
                    int s1 = (int)(7 * hz1);
                    int s2 = (int)(15 * hz1);
                    for (int j = 0; j < s1; j++)
                    {
                        mea_data[j].Set(0, 0);
                        mea_data[swp_data.Length - j - 1].Set(0, 0);
                    }
                    for (int j = s1; j < s2; j++)
                    {
                        double n = (double)(j - s1) / (s2 - s1);
                        double m = 0.5 * (1 + Math.Cos(Math.PI * (1 + n)));
                        mea_data[j].mul(m);
                        mea_data[swp_data.Length - j - 1].mul(m);
                    }
                }

                // Divide in complex domain
                for (int j = 0; j < swp_data.Length; j++)
                {
                    swp_data[j].idiv(mea_data[j]);
                    // Make a copy in mea_data, we'll need it later
                    mea_data[j] = swp_data[j];
                }

                // IFFT to get the impulse response
                Fourier.IFFT(swp_data.Length, swp_data);

                // Scan the imp to find maximum
                mx = 0;
                int mp = 0;
                for (int j = 0; j < sampleRate; j++)
                {
                    Complex si = swp_data[j];
                    double mg = Math.Abs(si.Magnitude);
                    if (mg > mx) { mx = mg; mp = j; }
                }
                // Look one sample backwards from max position
                // to find the likely "pulse peak" if within 30% of max
                peakpos = mp;
                if (mp>0 && swp_data[mp - 1].Magnitude / mx > 0.7)
                {
                    peakpos = mp - 1;
                }
            }

            // stderr.WriteLine("  {0} range {1}Hz to {2}Hz", isEstimatedSweepRange ? "Estimated sweep" : "Sweep", _fmin, _fmax);
            if (_fmaxSpecified && _fminSpecified)
            {
                HackSweep(swp_data, mea_data, peakpos, L, sampleRate);
            }
            else
            {
                Fourier.FFT(swp_data.Length, swp_data);
            }

            // Window the extremities of the whole response, finally?

            if (true)
            {
                // Make a copy in mea_data, we'll need it later
                for (int j = 0; j < swp_data.Length; j++)
                {
                    mea_data[j] = swp_data[j];
                }

                // Return FFT of impulse
                impulseFFT = mea_data;
            }

            // IFFT to get the impulse response
            Fourier.IFFT(swp_data.Length, swp_data);

            // Scan the imp to find maximum
            mx = 0;
            for (int j = 0; j < sampleRate; j++)
            {
                Complex si = swp_data[j];
                double mg = Math.Abs(si.Magnitude);
                if (mg > mx) mx = mg;
            }

            if (_noNorm)
            {
                mx = 1.0;
            }

            // Yield the first half (normalized) as result
            cs = new CallbackSource(1, sampleRate, delegate(long j)
            {
                if (j > (_returnAll ? L-1 : L / 2))
                {
                    return null;
                }
                Complex si = swp_data[j];
                Sample s = new Sample(si.Re / mx);
                return s;
            });
            cs.SampleRate = sampleRate;
            cs.NumChannels = 1;

            return cs;
        }

        static double mag(uint sampleRate, double f, double v)
        {
            // return Math.Log10(v);
            if (f > sampleRate / 2) f = sampleRate - f;
            double dbNow = 1 - (10 * Math.Log10(f));
            double gainNow = MathUtil.gain(dbNow);
            return v / gainNow;
        }
        static void get_stddev(uint sampleRate, int n1, int n2, Complex[] cdata, out double max, out double maxLog, out double stdev, out double avgLog)
        {
            int n3 = n2 - n1;
            max = 0;
            maxLog = 0;
            double tot = 0;
            double mean = 0;
            // Compute max magnitude, and totals
            for (int n = n1; n < n2; n++)
            {
                double v = cdata[n].Magnitude;
                if (v > max) max = v;
                double f = (double)n * sampleRate / cdata.Length;
                double val = mag(sampleRate, f, v);
                if (val > maxLog) maxLog = val;
                tot += val;
            }
            // Compute mean
            mean = tot / n3;
            // Compute deviation
            tot = 0;
            avgLog = 0;
            for (int n = n1; n < n2; n++)
            {
                double v = cdata[n].Magnitude;
                double f = (double)n * sampleRate / cdata.Length;
                double val = mag(sampleRate, f, v);
                tot += ((mean - val) * (mean - val));
                avgLog += val;
            }
            stdev = Math.Sqrt(tot / n3) / max;
            avgLog = avgLog / (n2 - n1);
        }

        static void HackSweep(Complex[] swp_data, Complex[] mea_data, int peakpos, int L, uint sampleRate)
        {
            int kmax = (int)(L / ((double)sampleRate / _fmax));
            int kmin = (int)(L / ((double)sampleRate / _fmin));

            // Construct FFT of a synthetic Dirac pulse at peakpos
            for (int j = 0; j < swp_data.Length; j++)
            {
                swp_data[j] = new Complex();
            }
            swp_data[peakpos] = new Complex(1, 0);
            Fourier.FFT(swp_data.Length, swp_data);
            // Calculate amplitude to scale the dirac pulse by
            double avg = 0;
            double avgI = 0;
            double gainfactor;
            for (int j = kmin; j < kmax; j++)
            {
                avg += (mea_data[j].Magnitude);   // sweep
                avgI += (swp_data[j].Magnitude);  // pulse
            }
            gainfactor = (avg/(kmax-kmin)) / (avgI/(kmax-kmin));

            // Window out the high frequencies beyond _fmax

            int _startSample = (int)((double)kmax * 0.98); //  k - 8000;
            int _endSample = (int)((double)kmax * 1.005); // k + 12000;
            double _startGain = 1;
            double _endGain = 0;
            int endData = _startSample;

            // Blackman-Harris window
            double c0 = 0.35875;
            double c1 = 0.48829;
            double c2 = 0.14128;
            double c3 = 0.01168;

            // Blackman
            c0 = 0.42; c1 = 0.5; c2 = 0.08; c3 = 0;

            for (int j = _startSample; j < _endSample; j++)
            {
                double frac = (double)(j - _startSample) / (double)(_endSample - _startSample);
                double phi = Math.PI * frac;
                double rcos = c0 - (c1 * Math.Cos(phi)) + (c2 * Math.Cos(2 * phi)) - (c3 * Math.Cos(3 * phi));
                double gain = _startGain + rcos * (_endGain - _startGain);

                Complex fill;

                fill = gainfactor * swp_data[j];
                swp_data[j] = (mea_data[j] * gain) + (fill * (1 - gain));

                fill = gainfactor * swp_data[L - j - 1];
                swp_data[L - j - 1] = (mea_data[L - j - 1] * gain) + (fill * (1 - gain));
            }
            for (int j = _endSample; j < L / 2; j++)
            {
                Complex fill;

                fill = gainfactor * swp_data[j];
                swp_data[j] = fill;

                fill = gainfactor * swp_data[L - j - 1];
                swp_data[L - j - 1] = fill;
            }

            // Window out the low frequencies
            // below _fmin
            _startSample = kmin / 2;
            _endSample = kmin;

            _startGain = 0;
            _endGain = 1;
            for (int j = 0; j < _startSample; j++)
            {
                Complex fill;

                fill = gainfactor * swp_data[j];
                swp_data[j] = fill;

                fill = gainfactor * swp_data[L - j - 1];
                swp_data[L - j - 1] = fill;
            }
            for (int j = _startSample; j < _endSample; j++)
            {
                double frac = (double)(j - _startSample) / (double)(_endSample - _startSample);
                double phi = Math.PI * frac;
                double rcos = c0 - (c1 * Math.Cos(phi)) + (c2 * Math.Cos(2 * phi)) - (c3 * Math.Cos(3 * phi));
                double gain = _startGain + rcos * (_endGain - _startGain);

                Complex fill;
                fill = gainfactor * swp_data[j];
                swp_data[j] = (mea_data[j] * gain) + (fill * (1 - gain));

                fill = gainfactor * swp_data[L - j - 1];
                swp_data[L - j - 1] = (mea_data[L - j - 1] * gain) + (fill * (1 - gain));
            }

            // Copy the real data back for everything else
            for (int j = _endSample; j < endData; j++)
            {
                swp_data[j] = mea_data[j];
                swp_data[L - j - 1] = mea_data[L - j - 1];
            }
        }


        static void Prep(string infileL, string infileR, string stereoImpulseFile, string outFile)
        {
            // Input files are complex
            // 0=mag, 1=phase/pi (so it looks OK in a wave editor!)
            // FFTs of the room impulse response

            // Take two half-FFT-of-impulse WAV files
            // Average them, into an array

            int n;
            SoundBuffer buff;
            WaveWriter wri;
            //          NoiseGenerator noise;
            FastConvolver conv;

            /*
            // Convolve noise with the in-room impulse
            noise = new NoiseGenerator(NoiseType.WHITE_FLAT, 2, (int)131072, stereoImpulse.SampleRate, 1.0);
            conv = new FastConvolver();
            conv.impulse = stereoImpulse;
            conv.Input = noise;

            wri = new WaveWriter("ImpulseResponse_InRoom.wav");
            wri.Input = conv;
            wri.Format = WaveFormat.IEEE_FLOAT;
            wri.BitsPerSample = 32;
            wri.Normalization = 0;
            wri.Run();
            wri.Close();
            wri = null;
            conv = null;
            */

            WaveReader rdrL = new WaveReader(infileL);
            buff = new SoundBuffer(rdrL);
            n = (int)buff.ReadAll();
            uint sampleRate = buff.SampleRate;
            uint nyquist = sampleRate / 2;

            double binw = (nyquist / (double)n);

            WaveReader rdrR = new WaveReader(infileR);

            IEnumerator<ISample> enumL = buff.Samples;
            IEnumerator<ISample> enumR = rdrR.Samples;

            // For easier processing and visualisation
            // read this in to an ERB-scale (not quite log-scale) array
            // then we can smooth by convolving with a single half-cosine.
            //

            int nn = (int)ERB.f2bin(nyquist, sampleRate) + 1;

            double[] muff = new double[nn];

            int prevbin = 0;
            int nbin = 0;
            double v = 0;
            int j = 0;
            while (true)
            {
                double f = (double)j * binw;    // equiv freq, Hz

                int bin = (int)ERB.f2bin(f, sampleRate); // the bin we drop this sample in
                if (bin > nn)
                {
                    // One of the channels has more, but we're overrun so stop now
                    break;
                }

                j++;
                bool more = false;
                more |= enumL.MoveNext();
                more |= enumR.MoveNext();
                if (!more)
                {
                    muff[prevbin] = v / nbin;
                    break;
                }

                v += enumL.Current[0];  // magnitude
                v += enumR.Current[0];  // magnitude
                nbin++;

                if (bin > prevbin)
                {
                    muff[prevbin] = v / nbin;
                    v = 0;
                    nbin = 0;
                    prevbin = bin;
                }
            }


            double[] smoo = ERB.smooth(muff, 38);

            // Pull out the freq response at ERB centers
            FilterProfile lfg = ERB.profile(smoo, sampleRate);

            // Write this response as a 'target' file
            /*
            FileStream fs = new FileStream("target_full.txt", FileMode.Create);
            StreamWriter sw = new StreamWriter(fs, Encoding.ASCII);
            foreach (FreqGain fg in lfg)
            {
                sw.WriteLine("{0} {1:f4}", Math.Round(fg.Freq), fg.Gain);
            }
            sw.Close();
            */
            /*
            fs = new FileStream("target_half.txt", FileMode.Create);
            sw = new StreamWriter(fs, Encoding.ASCII);
            foreach (FreqGain fg in lfg)
            {
                sw.WriteLine("{0} {1:f4}", Math.Round(fg.Freq), fg.Gain/2);
            }
            sw.Close();
            */


            // Create a filter to invert this response
            FilterProfile ifg = new FilterProfile();
            foreach (FreqGain fg in lfg)
            {
                ifg.Add(new FreqGain(fg.Freq, -fg.Gain));
            }

            ISoundObj filterImpulse = new FilterImpulse(0, ifg, FilterInterpolation.COSINE, sampleRate);
            filterImpulse.SampleRate = sampleRate;

            // Write the filter impulse to disk
            string sNoCorr = outFile + ".wav";
            wri = new WaveWriter(sNoCorr);
            wri.Input = filterImpulse; // invertor;
            wri.Format = WaveFormat.IEEE_FLOAT;
            wri.BitsPerSample = 32;
            wri.SampleRate = _sampleRate;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();
            _filterFiles.Add(sNoCorr);

            /*
            // Convolve noise with the NoCorrection filter
            noise = new NoiseGenerator(NoiseType.WHITE_FLAT, 2, (int)131072, stereoImpulse.SampleRate, 1.0);
            conv = new FastConvolver();
            conv.impulse = invertor;
            conv.Input = noise;

            wri = new WaveWriter("NoCorrection_Test.wav");
            wri.Input = conv;
            wri.Format = WaveFormat.IEEE_FLOAT;
            wri.BitsPerSample = 32;
            wri.SampleRate = _sampleRate;
            wri.Normalization = 0;
            wri.Run();
            wri.Close();
            wri = null;
            conv = null;
            */

            // Convolve this with the in-room impulse response

            WaveReader rea = new WaveReader(outFile + ".wav");
            conv = new FastConvolver();
            conv.impulse = rea;
            conv.Input = new WaveReader(stereoImpulseFile);
            wri = new WaveWriter(outFile + "_TestConvolution.wav");
            wri.Input = conv;
            wri.Format = WaveFormat.PCM;
            wri.Dither = DitherType.TRIANGULAR;
            wri.BitsPerSample = 16;
            wri.SampleRate = _sampleRate;
            wri.Normalization = -1;
            wri.Run();
            wri.Close();
            rea.Close();
            wri = null;
            conv = null;
        }


        #region DRC stuff: shelling out to DRC and assembling the results

        static bool DoDRC(string infileL, string infileR, string impDirectL, string impDirectR, int peakPosL, int peakPosR, string stereoImpulse, string stereoImpulseF)
        {
            // Call drc for each .drc file in the current directory.

            bool ok = true;
            FileInfo[] drcfiles = new DirectoryInfo(_thisFolder).GetFiles("*.drc");
            foreach (FileInfo drcfile in drcfiles)
            {
                // Run DRC for this drc file
                string outfile = Path.GetFileNameWithoutExtension(drcfile.Name);
                if (ExecDRC(drcfile.Name, null, outfile, infileL, infileR, peakPosL, peakPosR, stereoImpulse))
                {
                    _filterFiles.Add(outfile + ".wav");
                }
                else
                {
                    ok = false;
                }

                if (impDirectL != null && impDirectR != null)
                {
                    if (ExecDRC(drcfile.Name, null /*"flat.txt"*/, outfile + "_direct", impDirectL, impDirectR, peakPosL, peakPosR, stereoImpulse))
                    {
                        _filterFiles.Add(outfile + "_direct.wav");
                    }
                    else
                    {
                        ok = false;
                    }
                }
            }

            if (!ok)
            {
                stderr.WriteLine();
                stderr.WriteLine("Not all DRC ran successfully.");
            }
            return ok;
        }

        static bool ExecDRC(string drcfile, string target, string outfile, string infileL, string infileR, int peakPosL, int peakPosR, string stereoImpulseFile)
        {
            GC.Collect();
            bool ok = false;
            string args;
            FastConvolver conv;
            WaveWriter wri;
            if (!File.Exists(drcfile))
            {
                stderr.WriteLine();
                stderr.WriteLine("{0} not found.", drcfile);
                return ok;
            }

            string tmpL;
            string tmpR;
            tmpL = Path.GetFileNameWithoutExtension(infileL) + ".tmp";
            tmpR = Path.GetFileNameWithoutExtension(infileR) + ".tmp";
            _tempFiles.Add(tmpL);
            _tempFiles.Add(tmpR);

            stderr.WriteLine("Exec DRC for {0}, left channel", drcfile);
            stderr.WriteLine();
            ok = RunDRC(drcfile, infileL, target, tmpL, peakPosL, out args);

            if (ok)
            {
                stderr.WriteLine();
                stderr.WriteLine("Exec DRC for {0}, right channel", drcfile);
                stderr.WriteLine();
                ok = RunDRC(drcfile, infileR, target, tmpR, peakPosR, out args);
            }

            if (ok)
            {
                stderr.WriteLine();
                if (_noSkew)
                {
                    stderr.WriteLine("Creating stereo filter {0}", outfile + ".wav" );
                }
                else
                {
                    stderr.WriteLine("Creating stereo filter {0} (skew {1} samples)", outfile + ".wav", peakPosR - peakPosL);
                }
                ISoundObj stereoFilter = Splice(tmpL, peakPosL, tmpR, peakPosR, outfile + ".wav");
                stderr.WriteLine();

                // Convolve noise with the stereo filter
                /*
                NoiseGenerator noise = new NoiseGenerator(NoiseType.WHITE_FLAT, 2, (int)131072, stereoFilter.SampleRate, 1.0);
                conv = new FastConvolver();
                conv.impulse = stereoFilter;
                conv.Input = noise;

                wri = new WaveWriter(drcfile + "_Test.wav");
                wri.Input = conv;
                wri.Format = WaveFormat.IEEE_FLOAT;
                wri.BitsPerSample = 32;
                wri.SampleRate = _sampleRate;
                wri.Normalization = -1;
                wri.Run();
                wri.Close();
                wri = null;
                conv = null;
                noise = null;
                 * */

                // Convolve filter with the in-room impulse response
                WaveReader rea = new WaveReader(stereoImpulseFile);
                conv = new FastConvolver();
                conv.impulse = rea;
                conv.Input = stereoFilter;

                if (_pcm)
                {
                    _impulseFiles.Add("LCorrected_" + outfile + ".pcm: corrected test convolution, raw data (32-bit float), left channel");
                    _impulseFiles.Add("RCorrected_" + outfile + ".pcm: corrected test convolution, raw data (32-bit float), right channel");
                    WriteImpulsePCM(conv, "LCorrected_" + outfile + ".pcm", "RCorrected_" + outfile + ".pcm");
                }

                wri = new WaveWriter(outfile + "_TestConvolution.wav");
                wri.Input = conv;
                wri.Format = WaveFormat.PCM;
                wri.Dither = DitherType.TRIANGULAR;
                wri.BitsPerSample = 16;
                wri.SampleRate = _sampleRate;
                wri.Normalization = -1;
                wri.Run();
                wri.Close();
                rea.Close();
                wri = null;
                rea = null;
                conv = null;

                GC.Collect();
            }
            return ok;
        }

        static string GetDRCExe()
        {
            string exeName = "drc";
            if(File.Exists(Path.Combine(_thisFolder, "drc.exe")))
            {
                exeName = Path.Combine(_thisFolder, "drc.exe");
            }
            else if (File.Exists(Path.Combine(_thisFolder, "drc")))
            {
                exeName = Path.Combine(_thisFolder, "drc");
            }
            else if (DSPUtil.DSPUtil.IsMono() && File.Exists(Path.Combine("/usr/local/bin/", "drc")))
            {
                exeName = Path.Combine("/usr/local/bin/", "drc");
            }
            return exeName;
        }

        static string GetDRCArgs(string drcfile, string infile, string targetfile, string outfile, int peakPos)
        {
            string args = escapeArg(drcfile);
            args += " --BCInFile " + escapeArg(infile);
            args += " --BCInFileType F";                // input is floating-point
            args += " --BCImpulseCenter " + peakPos;    // we know the center of the impulse
            args += " --BCImpulseCenterMode M";         // M = specified, not auto
            args += " --MCOutFile " + escapeArg(outfile);          // this is our output file (after mic compensation)
            args += " --MCOutFileType F";               // output floating-point
            if (!_noOverrideDRC)
            {
                args += " --MCOutWindow " + _filterLen;     // length of output file
                args += " --MCFilterLen " + _filterLen;     // length of output file
                args += " --MCFilterType L";                // linear-phase mic compensation filter
            }
            if (targetfile != null)
            {
                args += " --PSPointsFile " + escapeArg(targetfile);    // target response file
                args += " --PSMagType D";                   // target magnitudes are dB
            }
            if (!_noOverrideDRC)
            {
                args += " --PSFilterType L";                // linear-phase target filter
                args += " --PSOutWindow " + _filterLen;     // size of the output file
                args += " --MSOutWindow " + _filterLen;     // in case this is used
            }
            return args;
        }

        static string escapeArg(string arg)
        {
            if (arg.Contains(" "))
            {
                arg = "\"" + arg + "\"";
            }
            return arg;
        }

        static bool RunDRC(string drcfile, string infile, string targetfile, string outfile, int peakPos, out string drcArgs)
        {
            string exeName = escapeArg(GetDRCExe());

            drcArgs = GetDRCArgs(drcfile, infile, targetfile, outfile, peakPos);
            stderr.WriteLine(drcArgs);

            System.Diagnostics.Process drcProcess = new System.Diagnostics.Process();
            System.Diagnostics.ProcessStartInfo drcInfo = new System.Diagnostics.ProcessStartInfo();
            drcInfo.Arguments = drcArgs;
            drcInfo.FileName = exeName;
            drcInfo.ErrorDialog = false;
            drcInfo.CreateNoWindow = true;
            drcInfo.UseShellExecute = false;
            drcInfo.RedirectStandardOutput = true;

            drcProcess.StartInfo = drcInfo;

            bool ok = false;
            try
            {
                stderr.WriteLine();
                drcProcess.Start();
                while (!drcProcess.HasExited)
                {
                    while (!drcProcess.StandardOutput.EndOfStream && !drcProcess.HasExited)
                    {
                        stderr.WriteLine(drcProcess.StandardOutput.ReadLine());
                    }
                    drcProcess.WaitForExit(1000);
                }
                stderr.Write(drcProcess.StandardOutput.ReadToEnd());
                stderr.WriteLine();
                ok = true;
            }
            catch (Exception e)
            {
                stderr.WriteLine("DRC failed to run: " + e.Message);
            }
            if (drcProcess.ExitCode != 0)
            {
                stderr.WriteLine("DRC failed.  Command: " + exeName + " " + drcArgs);
                ok = false;
            }
            return ok;
        }

        static ISoundObj Splice(string infileL, int peakPosL, string infileR, int peakPosR, string outfile)
        {
            //int tmp1;
            //bool tmp2;
            WaveReader reader1 = new WaveReader(infileL, WaveFormat.IEEE_FLOAT, 32, 1);
//            reader1.Skip(peakPosL, out tmp1, out tmp2);
            SoundBuffer buff1 = new SoundBuffer(reader1);
            buff1.ReadAll();

            // Normalize loudness for each channel
            double g = Loudness.WeightedVolume(buff1);
            buff1.ApplyGain(1/g);

            WaveReader reader2 = new WaveReader(infileR, WaveFormat.IEEE_FLOAT, 32, 1);
//            reader2.Skip(peakPosR, out tmp1, out tmp2);
            SoundBuffer buff2 = new SoundBuffer(reader2);
            buff2.ReadAll();

            g = Loudness.WeightedVolume(buff2);
            buff2.ApplyGain(1/g);

            ChannelSplicer splicer = new ChannelSplicer();
            splicer.Add(buff1);
            splicer.Add(buff2);

            // General-purpose:
            // Find the extremities of the DRC impulse,
            // window asymmetrically.
            //
            // But, since we specifically used linear-phase filters on target and mic,
            // we know that the impulse is centered.
            // We want an impulse length (_filterLen)
            // so window the 2048 samples at each end (which are pretty low to begin with - less than -100dB)

            ISoundObj output;
            int nCount = (int)(buff1.Count / 2);
            if (nCount > _filterLen/2)
            {
                BlackmanHarris bhw = new BlackmanHarris(nCount, 2048, (nCount / 2) - 2048);
                bhw.Input = splicer;
                SampleBuffer sb = new SampleBuffer(bhw);
                output = sb.Subset(_filterLen/2, _filterLen);
            }
            else
            {
                output = splicer;
            }

            ISoundObj result = output;
            if (!_noSkew)
            {
                // Apply skew to compensate for time alignment in the impulse responses
                Skewer skewer = new Skewer(true);
                skewer.Input = output;
                skewer.Skew = peakPosL - peakPosR;
                result = skewer;
            }

            WaveWriter writer = new WaveWriter(outfile);
            writer.Input = result;
            writer.Dither = DitherType.NONE;
            writer.Format = WaveFormat.IEEE_FLOAT;
            writer.BitsPerSample = 32;
            writer.SampleRate = _sampleRate;
            writer.NumChannels = splicer.NumChannels;
            if (Double.IsNaN(_gain))
            {
                writer.Normalization = 0;
            }
            else
            {
                writer.Gain = MathUtil.gain(_gain);
            }
            writer.Run();
            writer.Close();
            reader1.Close();
            reader2.Close();

            return result;
        }

        #endregion


        /// <summary>
        /// Output pipe
        /// </summary>
        static TextWriter stderr
        {
            get
            {
                return System.Console.Out;
            }
        }


        #region Usage and help

        /// <summary>
        /// Show app name and version
        /// </summary>
        /// <returns>false if expired</returns>
        static bool DisplayInfo()
        {
            bool ok = true;
            System.Console.Error.WriteLine("ImpulsePrep v{0}", DSPUtil.DSPUtil.VERSION);
            System.Console.Error.WriteLine("http://inguzaudio.com/Tools/ImpulsePrep/");
            /*
            if (DSPUtil.DSPUtil.EXPIRY != null)
            {
                if (DateTime.Now.CompareTo(DSPUtil.DSPUtil.EXPIRY) >= 0)
                {
                    System.Console.Error.WriteLine("**** THIS EVALUATION VERSION EXPIRED {0}", DSPUtil.DSPUtil.EXPIRY);
                    ok = false;
                }
                else
                {
                    System.Console.Error.WriteLine("This evaluation version will expire {0}", DSPUtil.DSPUtil.EXPIRY);
                }
            }
            */
            stderr.WriteLine();
            return ok;
        }


        /// <summary>
        /// Show the help
        /// </summary>
        /// <param name="longUsage"></param>
        static void DisplayUsage(bool longUsage)
        {
            stderr.WriteLine("Usage: ImpulsePrep /L filename /R filename [/dbl] [/split] [/nodrc]");
            stderr.WriteLine("          [/direct] [/nosub] [/copy] [/length n] [/noover] [/pcm] [/refch n]");
            stderr.WriteLine("");
            stderr.WriteLine("  /L : WAV file with room recording for left channel");
            stderr.WriteLine("  /R : WAV file with room recording for right channel");
            stderr.WriteLine("");
            stderr.WriteLine("Each recording file should be stereo;");
            stderr.WriteLine("- one channel with a direct line recording of the sweep,");
            stderr.WriteLine("- second channel with microphone-recorded response.");
            stderr.WriteLine("");
            stderr.WriteLine("The left and right recordings are processed to calculate the impulse response");
            stderr.WriteLine("of your amplifier, loudspeakers and room.");
            stderr.WriteLine("");
            stderr.WriteLine("The left and right impulse responses are usually then automatically processed");
            stderr.WriteLine("to create room-correction filters using Denis Sbragion's DRC");
            stderr.WriteLine("(http://drc-fir.sourceforge.net/).");
            stderr.WriteLine("");
            if (!longUsage)
            {
                stderr.WriteLine("For detailed help: ImpulsePrep /?");
            }
            else
            {
                stderr.WriteLine("");
                stderr.WriteLine("Optional parameters:");

                stderr.WriteLine("");
                stderr.WriteLine("/DBL");
                stderr.WriteLine("    With this flag, the left and right impulses are also saved in .DBL format");
                stderr.WriteLine("    ready for processing with Acourate (http://www.acourate.com/).");

                stderr.WriteLine("");
                stderr.WriteLine("/SPLIT");
                stderr.WriteLine("    With this flag, the sweep and recording are separated and saved as PCM,");
                stderr.WriteLine("    suitable for processing with DRC's 'lsconv' utility.");

                stderr.WriteLine("");
                stderr.WriteLine("/NODRC");
                stderr.WriteLine("    This flag disables DRC for filter creation.");

                stderr.WriteLine("");
                stderr.WriteLine("/DIRECT");
                stderr.WriteLine("    Additionally creates a filtered version of the impulse response,");
                stderr.WriteLine("    approximately corresponding to the direct sound (progressively filtering to");
                stderr.WriteLine("    remove the reverberant field).  This 'direct' response is used to make a");
                stderr.WriteLine("    second, alternate set of room correction filters.");

                stderr.WriteLine("");
                stderr.WriteLine("/NOSUB");
                stderr.WriteLine("    By default, ImpulsePrep applies a gentle subsonic filter (below 15Hz) to");
                stderr.WriteLine("    the measured impulse response, which avoids very low frequency noise.");
                stderr.WriteLine("    With the /NOSUB flag, subsonic filtering is skipped.");

                stderr.WriteLine("");
                stderr.WriteLine("/REFCH");
                stderr.WriteLine("    Specifies that channel <N> is the reference channel (the original sweep)");
                stderr.WriteLine("    (0=left, 1=right).  Only needed when automatic detection of the reference");
                stderr.WriteLine("    channel doesn't work and the resulting impulses look like pure noise.");

                stderr.WriteLine("");
                stderr.WriteLine("/GAIN");
                stderr.WriteLine("    Applies gain (dB) to the correction filters.");

                stderr.WriteLine("");
                stderr.WriteLine("");
                stderr.WriteLine("Optional parameters when DRC is used:");

                stderr.WriteLine("");
                stderr.WriteLine("/COPY");
                stderr.WriteLine("    With this flag, correction filters created by DRC are copied into the");
                stderr.WriteLine("    Impulses folder for InguzDSP:");
                stderr.WriteLine("    " + _impulsesFolder);

                stderr.WriteLine("");
                stderr.WriteLine("/LENGTH");
                stderr.WriteLine("    Specifies a length for the correction filters produced by DRC");
                stderr.WriteLine("    (default: 32768)");

                stderr.WriteLine("");
                stderr.WriteLine("/NOOVER");
                stderr.WriteLine("    With this flag, your .drc file's settings for MCOutWindow, MCFilterLen,");
                stderr.WriteLine("    MCFilterType, PSFilterType, PSOutWindow and MSOutWindow are not overridden");
                stderr.WriteLine("    when ImpulsePrep calls DRC.");

                stderr.WriteLine("");
                stderr.WriteLine("/NOSKEW");
                stderr.WriteLine("    This flag disables interchannel time-alignment of the correction filters.");

                stderr.WriteLine("");
                stderr.WriteLine("/PCM");
                stderr.WriteLine("    With this flag, left and right impulses and corrected (test convolution)");
                stderr.WriteLine("    impulses are also saved in .PCM format ready for graphing with the");
                stderr.WriteLine("    Octave scripts supplied with DRC");
                stderr.WriteLine("    (http://drc-fir.sourceforge.net/doc/drc.html#htoc207).");

                stderr.WriteLine("");
            }
            stderr.WriteLine("");
        }

        #endregion

    }
}

This fork of [gnss-sdr](https://github.com/oscimp/gnss-sdr) aims at 
providing spoofing detection capability by analyzing the direction of 
arrival of the signals transmitted from each GPS satellite transmitting 
in the L1 band. It is assumed that two antennas are connected to the
two inputs of a dual channel coherent SDR receiver -- tests were completed with
the Ettus Research B210 -- separated by half a wavelength (10 cm at L1).

The original gnss-sdr installation documentation is found in [README.original](README.original).

This software was tested with an Ettus Research B210 dual-input SDR platform and
with the File Source.

The Signal_Source philosophy is probably broken by including the spoofing detection
processing in [spoofing_detection](src/algorithms/signal_source/libs/spoofing_detection.cc). 
This solves the issue of multiple antenna-inputs and single output.

Activating the Spoofing Detection flag is achieved by adding the
``SignalSource.spoofing_protection=1`` flag in the configuration file.

Most basic compilation:
```shell
cd build
cmake ../
make -j4
./src/main/gnss-sdr -c ../File-GNSS-SDR-receiver.conf
```
assuming the ``File-GNSS-SDR-receiver.conf`` has been updated to point to an existing
pair of files recorded e.g. from a B210 as a spoofing signal was being emitted. The
File Source format is to provide directory + beginning of the file name and the
extension '_1.bin' and '_2.bin' will be added when loading the files. In this example,
the files are
```shell
$ ls -l /t/7_m35dBm/
-rw-r--r-- 1 xxx xxx 1964160000 Jun  5 17:23 7_m35dBm_1.bin
-rw-r--r-- 1 xxx xxx 1964160000 Jun  5 17:24 7_m35dBm_2.bin
```

Running the spoofing detection mechanism from gnss-sdr on these files will display
```shell
10:     meanarg=0.6722  meanabs=8.122   stdargres_=0.00061      weightabs=8.09,weightarg=0.68 /!\
```
with stdargres_=0.00061 meaning spoofing is occuring (the standard deviation on the angle
of arrival is too low to be compatible with a genuine constellation). The threshold detection
indicating spoofing was triggered with the /!\ sign at the end of the line.

On a genuine constellation,
```shell
13:     meanarg=-0.2871 meanabs=7.054   stdargres_=3.39117      weightabs=0.00,weightarg=0.00
```
indicates, through its large stdargres_ value, that no spoofing is occuring. Under such conditions,
decoding is performed as would be done with a classical gnss analysis sequence
```shell
8:      meanarg=-0.5817 meanabs=6.539   stdargres_=3.87295      weightabs=0.00,weightarg=0.00
New GPS NAV message received in channel 0: subframe 3 from satellite GPS PRN 01 (Block IIF)
New GPS NAV message received in channel 17: subframe 3 from satellite GPS PRN 22 (Block IIR)
New GPS NAV message received in channel 19: subframe 3 from satellite GPS PRN 14 (Block IIR)
New GPS NAV message received in channel 10: subframe 3 from satellite GPS PRN 17 (Block IIR-M)
New GPS NAV message received in channel 12: subframe 3 from satellite GPS PRN 32 (Block IIF)
New GPS NAV message received in channel 9: subframe 3 from satellite GPS PRN 28 (Block IIR)
First position fix at 2019-Nov-26 07:58:00.120000 UTC is Lat = 47.2534 [deg], Long = 5.99282 [deg], Height= 540.828 [m]
7:      meanarg=-0.3519 meanabs=7.124   stdargres_=3.90482      weightabs=0.00,weightarg=0.00
Position at 2019-Nov-26 07:58:00.500000 UTC using 4 observations is Lat = 47.254877320 [deg], Long = 5.994379292 [deg], Height = 885.775 [m]
Velocity: East: -0.045 [m/s], North: 0.182 [m/s], Up = 0.434 [m/s]
```

A graphical representation of this result is shown of the figure below

<img src="sortie_phase_zeropm2suitesuite.png">

where a genuine record was collected (top=Doppler shift of the detected satellite, middle=phase
between antennas, bottom=position) initially (left), then with a spoofing signal generating
a location West of France in Britanny (Brest) for 6 minutes, before returning to the genuine 
signal, and finally genuine signal decoding (right) of the correct location in Besancon (East of
France).

See 
[1] and [2] for an explanation on the analaysis of the standard deviation of the phase between
antennas.

[1] J.-M. Friedt, W. Feng
Anti-leurrage et anti-brouillage de GPS par réseau d'antennes, MISC 110 (2020) [in French]

[2] J.-M. Friedt, W. Feng, G. Goavec-Merou, F. Meyer, GPS spoofing implementation by Software Defined 
Radio and computationally efficient GPS spoofing detection and cancellation (submitted, 2020)


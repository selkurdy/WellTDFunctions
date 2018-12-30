##  WellTDf

##  CONCEPT

>  *__welltdf__* reads in data either from a 2 column depth time files or a las like flat file with the depth and dt, i.e. slowness values. These are then converted to depth one way time pairs and a quadratic function is fitted. Also, the instantaneous, in case of sonic input, or interval, in case of depth time tables input, velocities are computed and a linear function `V0,K` is also computed. 

>  These coefficients per well are then accummulated in another text file to be mapped and used for depth conversion.


##  Command Line

```
python welltdf.py -h
usage: welltdf.py [-h] [--dtcols DTCOLS DTCOLS] [--headerlines HEADERLINES]
                  [--slowness] [--tometric] [--nulval NULVAL]
                  [--dtshift DTSHIFT DTSHIFT]
                  [--dtmultiplier DTMULTIPLIER DTMULTIPLIER]
                  [--topsfilename TOPSFILENAME] [--topcols TOPCOLS TOPCOLS]
                  [--topsin {time,depth}]
                  [--topsmultshift TOPSMULTSHIFT TOPSMULTSHIFT]
                  [--wellname WELLNAME] [--xycoords XYCOORDS XYCOORDS]
                  [--zrange ZRANGE ZRANGE ZRANGE] [--listdatain] [--listdt]
                  [--hideplot]
                  datafilename

Fit Linear with Depth and Quadratic functions to sonic or DTR. Dec 30,2018

positional arguments:
  datafilename          Depth Time 1Way file with sonic, density,shear

optional arguments:
  -h, --help            show this help message and exit
  --dtcols DTCOLS DTCOLS
                        depth T1W columns. dfv=0 1
  --headerlines HEADERLINES
                        headerlines to skip. dfv=0
  --slowness            Input is LAS, i.e. depth slowness. dfv= depth time
                        pairs
  --tometric            convert slowness to velocities in m/s.default= keep as
                        ft/s
  --nulval NULVAL       Null value in sonic log. dfv= -999.25000
  --dtshift DTSHIFT DTSHIFT
                        depth shift time shift in ms, dfv = 0 0
  --dtmultiplier DTMULTIPLIER DTMULTIPLIER
                        depth and time multipliers
  --topsfilename TOPSFILENAME
                        file name with 2 columns, depth and tops
  --topcols TOPCOLS TOPCOLS
                        columns of depth and top,dfv= 0 1
  --topsin {time,depth}
                        tops are in depth or t2w. dfv= depth
  --topsmultshift TOPSMULTSHIFT TOPSMULTSHIFT
                        Tops mutliplier and shift
  --wellname WELLNAME   Well name column, dfv=WELL
  --xycoords XYCOORDS XYCOORDS
                        X Y coordinates of well.dfv=664055.00 ,2889245.00
  --zrange ZRANGE ZRANGE ZRANGE
                        Depth range to fit functions and plot. Expected 0/1
                        zmin zmax values. dfv= 0 0 0
  --listdatain          List input data
  --listdt              List d t data
  --hideplot            Only save to pdf. default =show and save

```  


>  `welldtf` expects a file name which has to be given on the command line. The file is expected to be a flat columnar file, i.e. __no csv__   
>  `--headerlines` # of header lines to skip to reach the start of the data.   
>  `--dtcols` 2 integer values are expected to represent the depth and the time column numbers. The default is 0 1   
>  `--slowness`  use this option to indicate that the input time is sonic dt. the default is one way time in seconds.  
>  `--tometric` is used only with `--slowness`  option to indicate that you are working in meters/second for velocities.  
>  `--dtmultiplier` 2 values are expected, the first is for depth and the second is for time. This option is useful to convert time from two way to oneway and from ms to sec, e.g. 1.0 0.0005  
>  `--zrange` use this option to select a subset of the data by depth range.  



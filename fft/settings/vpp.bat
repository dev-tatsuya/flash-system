set path_1=%path%
set MSDevDir=c:\Visual_Studio\Vc98
set LIB=%MSDevDir%\lib;
path c:\Visual_Studio\Vc98\bin;c:\Visual_Studio\Vc98\include;c:\Visual_Studio\Common\MsDev98\bin
CL /Ic:\Visual_Studio\SLSWIN /I%MSDevDir%\INCLUDE %1 %2 %3 %4 %5 %6 %7 %8 %9 /link user32.lib gdi32.lib
set LIB=
path ;
path %path_1%


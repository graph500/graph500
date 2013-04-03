rem This does not work because VS2012 is not fully C++11 compliant yet.
setlocal
set CFLAGS=/DR123_USE_CXX11=1 /I..\include /W3 /Ox /EHs /nologo /favor:INTEL64
set VCBAT="c:\Program Files (x86)\Microsoft Visual Studio 11.0\vc\vcvarsall.bat"
call BUILDVC.bat %1
endlocal

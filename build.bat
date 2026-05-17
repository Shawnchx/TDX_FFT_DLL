@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvars64.bat"
"C:\Program Files\Microsoft Visual Studio\18\Insiders\MSBuild\Current\Bin\MSBuild.exe" D:\TDX_FFT_DLL\TDX_FFT_DLL.sln /p:Configuration=Release /p:Platform=x64 /p:PlatformToolset=v143 /t:Build /verbosity:normal

@echo off
echo ============================================
echo  Building TDX_FFT_DLL v5 (phase: t=0 at newest bar)
echo ============================================
call "C:\Program Files\Microsoft Visual Studio\18\Insiders\VC\Auxiliary\Build\vcvars64.bat"
"C:\Program Files\Microsoft Visual Studio\18\Insiders\MSBuild\Current\Bin\MSBuild.exe" D:\TDX_FFT_DLL\v5\TDX_FFT_DLL_v5.vcxproj /p:Configuration=Release /p:Platform=x64 /p:PlatformToolset=v143 /t:Build /verbosity:normal
echo.
echo Output: D:\TDX_FFT_DLL\v5\bin\x64\Release\TDX_FFT_DLL_v5.dll
echo Deploy to TDX as: T0002\dlls\TDX_FFT_DLL_v5.dll  (bind as slot 2)

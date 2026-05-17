@echo off
setlocal enabledelayedexpansion

set "VCTOOLS=C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build"

if not exist "%VCTOOLS%\vcvarsall.bat" (
    echo [错误] 未找到 VS2022，尝试查找 VS2019...
    set "VCTOOLS=C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Auxiliary\Build"
    if not exist "!VCTOOLS!\vcvarsall.bat" (
        echo [错误] 未找到 Visual Studio 编译环境，请确认已安装。
        echo 如果已安装，请手动打开"开发者命令提示符"后运行 msbuild。
        exit /b 1
    )
)

echo [编译] 配置环境: x64
call "%VCTOOLS%\vcvarsall.bat" x64

echo [编译] 正在编译 v6 Release x64...
msbuild TDX_FFT_DLL_v6.vcxproj /p:Configuration=Release /p:Platform=x64

if %errorlevel% equ 0 (
    echo [成功] 编译完成: bin\x64\Release\TDX_FFT_DLL_v6.dll
) else (
    echo [失败] 编译出错，请检查代码。
)
pause

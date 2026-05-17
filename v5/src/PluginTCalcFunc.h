// PluginTCalcFunc.h
// 通达信DLL插件接口规范 - 原版官方头文件（保持完全一致）

#ifndef __PLUGIN_TCALC_FUNC
#define __PLUGIN_TCALC_FUNC
#pragma pack(push,1) 

// 函数原型: (数据个数, 输出数组, 输入A, 输入B, 输入C)
typedef void(*pPluginFUNC)(int, float*, float*, float*, float*);

typedef struct tagPluginTCalcFuncInfo
{
    unsigned short  nFuncMark;   // 函数编号
    pPluginFUNC     pCallFunc;   // 函数地址
} PluginTCalcFuncInfo;

typedef BOOL(*pRegisterPluginFUNC)(PluginTCalcFuncInfo**);

#pragma pack(pop)
#endif

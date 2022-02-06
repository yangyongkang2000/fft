# 一维快速傅立叶变换
## 算法介绍
> 二次幂采用CooleyTukey算法，非二次幂采用Bluestein算法，针对小范围数据使用DFT。Bluestein算法比较Naive,虽然保证了O(N*Log(N))时间复杂度，但使用了额外空间，并且精度较低，运行效率远不如二次幂高效。

## 简要介绍
+ C++ 模版泛型编程 仅需要一个头文件
+  “API简单易用高效且API易扩展,你只要会用STL,就会用”。
+ 需要支持C++20以上的编译器（因为使用了\<bit> 头文件）


## 基础API介绍

```C++
template<bool parallel,bool Inverse,bool R,typename InputIt>
inline auto CooleyTukey(InputIt beg,InputIt end,const unsigned int N=0) noexcept
```

+ 这是专门针对2次幂的fft的裸函数


+ parallel代表是否并行，Inverse 代表是否是逆变换。R代表是否需要调整数组。InputIt为一维数组迭代器类型。N代表buffer位置，默认为0，主要用来控制不同的buffer位置,进行并行的时候可以避免数据竞争。


> 我来展示一下代码示例

```C++
#include <iostream>
#include<complex>
#include"fourier_transform.hpp"
int main(int argc, const char * argv[]) {
    // insert code here...
    std::complex<double> p[8] {1,2,3,4,5,6,7,8};
    yyk_fft::CooleyTukey<false, false, true>(p, p+8);
    for(auto &z:p) std::cout<<z<<" ";
    std::cout<<'\n';
    return 0;
}
```

> 当然你不一定需要标准库的complex，你也可以用自己造的complex，都没问题，设计的API对complex类的要求非常少。CooleyTukey不需要开辟额外空间，对输入迭代器范围进行原位傅立叶变换。



+ 当然还有专门针对非二次幂的Bluestein算法
```C++
template<bool parallel,bool Inverse,typename InputIt>
inline InputIt Bluestein(InputIt beg,InputIt end) noexcept
```
+ 也有O(N^2)的DFT
```C++
template<bool Inverse,typename InputIt>
inline InputIt dft(InputIt beg,InputIt end,const unsigned int n=0) noexcept
```


### 以上这些API并不建议直接裸用，因为计算fft,有封装好的API去直接调用,但它们也很重要，是原始的轮子或者说原始的物料，用来造新轮子的。




## 封装好后的API

### 变换结果没有归一化


+ 快速傅立叶变换
```C++
template<bool parallel=false,typename InputIt>
inline InputIt fast_fourier_transform(InputIt beg,InputIt end) noexcept
```

+ 使用起来超简单。示例
```C++
#include <iostream>
#include<complex>
#include<numeric>
#include<vector>
#include"fourier_transform.hpp"
int main(int argc, const char * argv[]) {
    // insert code here...
    std::vector<std::complex<double>> p(100);
    std::iota(p.begin(),p.end(), 1);
    yyk_fft::fast_fourier_transform(p.begin(),p.end());
    for(auto &z:p) std::cout<<z<<'\n';
    return 0;
}
```

+ 快速傅立叶逆变换

```C++
template<bool parallel=false,typename InputIt>
inline InputIt inverse_fast_fourier_transform(InputIt beg,InputIt end) noexcept
```

+ 使用起来超简单。示例
```C++
#include <iostream>
#include<complex>
#include<numeric>
#include<vector>
#include"fourier_transform.hpp"
int main(int argc, const char * argv[]) {
    // insert code here...
    std::vector<std::complex<double>> p(100);
    std::iota(p.begin(),p.end(), 1);
    yyk_fft::inverse_fast_fourier_transform(p.begin(),p.end());
    for(auto &z:p) std::cout<<z<<'\n';
    return 0;
}
```

## Wolfram语言接口

> 利用WolframLibraryLink技术，很容易将API编译成动态库给Wolfram使用

+ 以下是源文件，我在MacOS下编译成yyk_fft.dylib
```C++
#define CPP_Complex
#include <stdio.h>
#include<complex>
#include<algorithm>
#include<cmath>
#include<WolframLibrary.h>
#include"fourier_transform.hpp"
EXTERN_C DLLEXPORT  int fft(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) noexcept
{
    auto tensor=MArgument_getMTensor(Args[0]);
    auto b=MArgument_getBoolean(Args[1]);
    auto data=libData->MTensor_getComplexData(tensor);
    auto len=libData->MTensor_getFlattenedLength(tensor);
    if (b)
        yyk_fft::fast_fourier_transform<true>(data, data+len);
    else yyk_fft::fast_fourier_transform(data, data+len);
    MArgument_setMTensor(Res, tensor);
    return LIBRARY_NO_ERROR;
}
EXTERN_C DLLEXPORT  int ifft(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res) noexcept
{
    auto tensor=MArgument_getMTensor(Args[0]);
    auto b=MArgument_getBoolean(Args[1]);
    auto data=libData->MTensor_getComplexData(tensor);
    auto len=libData->MTensor_getFlattenedLength(tensor);
    if (b)
        yyk_fft::inverse_fast_fourier_transform<true>(data, data+len);
    else yyk_fft::inverse_fast_fourier_transform(data, data+len);
    MArgument_setMTensor(Res, tensor);
    return LIBRARY_NO_ERROR;
}
```
> 然后将动态库导入
```mathematica
g = LibraryFunctionLoad[
  "yyk_fft.dylib", 
  "ifft", {{Complex, 1}, True | False}, {Complex, 1}];
fft[l_] := f[l, False];
pfft[l_] := f[l, True];
ifft[l_] := g[l, False];
pifft[l_] := g[l, True];
```
> 虽然没有直接写二维傅立叶变换API，但稍微改一下，通过一维傅立叶变换间接操作就可以进行二维傅立叶变换

```mathematica
fft2[l_] := Transpose[fft /@ Transpose[fft /@ l]]~Divide~Length[l]
```

## 写到最后

+ 一些缺陷或者Bug
> 并行的缺陷

> 内存的缺陷

> 算法的缺陷，针对于非二次幂的算法。


> 功能的缺陷，比如没有二维傅立叶变换，虽然可以很容易在一维傅立叶变换基础上开发，但是我并没有，如果读者感兴趣，我会考虑，否则留给大家吧。

> 计算精度和计算结果是否正确，留给读者检验

> 还有很多... 希望读者提出改进建议和指出不足和缺陷以及Bug

+ 性能问题

> 性能是用C++提高了C抽象的同时保持住C的高效，还是抽象之后，性能的搞笑? 这些疑问留给读者评价。

> 优化全靠编译器

+ 奔跑吧,C++
```shell
 g++ -Ofast -march=native ...
```
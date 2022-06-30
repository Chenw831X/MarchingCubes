>《科学计算可视化》上机作业
>
>姓名：陈威
>
>学号：12121085

# Marching cubes

本项目实现了**Marching cubes**算法，并使用**并行**编程库OpenMP加速；使用[libigl](https://libigl.github.io/)导出OBJ格式文件。

* 顶点和边的索引使用:

<img src="/home/cw/.config/Typora/typora-user-images/image-20220629230706651.png" alt="image-20220629230706651" style="zoom:50%;" />

* lookup table参考[Polygonising a scalar field](http://lemur.cmp.uea.ac.uk/Research/ivis/backup/PhD/Ronan%20iViS%20Stuff%20%28Website%29/Polygonising%20a%20scalar%20field.pdf)

### 项目结构

```C++
.
├── Utils.hpp  // 读取二进制raw格式的数据，线性插值计算交点
├── MarchingCubes.hpp // Marching cubes, lookup-table, 导出OBJ格式文件
├── MarchingCubes.cpp
```

### 项目编译

#### 环境

* `Ubuntu 22.04`, `Intel i7-11700 2.5GHz`
* `cmake 3.22`, `g++ 11.2.0`

#### 依赖

1. 项目依赖于并行编程库OpenMP
2. 此外，项目还用到了：

   * `Eigen`：向量运算

   * `libigl`：导出OBJ格式文件

   * `spdlog`：用于日志输出


#### 编译

```shell
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```



### 结果展示

| iso value (data[i] $\in [0,4964]$) | time cost |
| ---------------------------------- | --------- |
| 1000                               | 1h27min   |
| 1861                               | 112s      |

#### 2. iso value = 1000

<img src="/home/cw/.config/Typora/typora-user-images/image-20220630091341117.png" alt="image-20220630091341117" style="zoom:50%;" />



#### 1. iso value = 1861

<img src="/home/cw/.config/Typora/typora-user-images/image-20220629224814199.png" alt="image-20220629224814199" style="zoom:50%;" />
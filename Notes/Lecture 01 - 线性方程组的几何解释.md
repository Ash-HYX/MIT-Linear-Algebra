---
title: "MIT 18.06 Lecture 01: The Geometry of Linear Equations"
course: "MIT 18.06 Linear Algebra (Gilbert Strang)"
book: "Introduction to Linear Algebra, Gilbert Strang"
tags: [线性代数, 18.06, MIT, 线性方程组, 线性组合, 行图像, 列图像, 矩阵形式, 奇异矩阵]
---

# Lecture 01 — 线性方程组的几何解释

> 本讲主线：线性代数的**根本问题**是求解线性方程组。同一个方程组有 **三种等价视角**——行图像、列图像、矩阵形式。其中**列图像（column picture）最重要**：求解 $A\mathbf{x}=\mathbf{b}$ 就是「用 $A$ 的**各列**线性组合出 $\mathbf{b}$」，解 $\mathbf{x}$ 就是那组组合系数。

本节只讨论「正常情形」：**$n$ 个方程、$n$ 个未知数**（方程个数 = 未知数个数，方形方程组）。

---

## 1. 基本问题与三种视角

- 问题：求解线性方程组（system of linear equations）。
- 三种视角（同一件事的三种说法）：
  1. **行图像（row picture）**：一次看一个方程。$n=2$ 时是两条直线相交，$n=3$ 时是三个平面交于一点。
  2. **列图像（column picture）**：一次看一列。把 $A$ 的列向量做**线性组合**去凑出右端向量 $\mathbf{b}$。
  3. **矩阵形式（matrix form）**：整体写成 $A\mathbf{x}=\mathbf{b}$，是紧凑的记号与运算载体。

> **为什么列图像是关键**：行图像的几何直观在 $n\ge 3$ 之后迅速失效（三维画三个平面已经很勉强，四维以上无法可视化）；而「列的线性组合」这一想法可以在任意维数下继续使用，是整个课程后续理论（列空间、秩、消元法）的载体。

---

## 2. 矩阵形式

$$
A\mathbf{x}=\mathbf{b}
$$

- $A$：**系数矩阵（coefficient matrix）**，由方程组中未知数的系数按行排列而成的矩形数组。
- $\mathbf{x}$：未知向量（unknown vector）。
- $\mathbf{b}$：右端项 / 右端向量（right-hand side vector）。

**维度规则**：$A$ 的行数 = 方程个数，$A$ 的列数 = 未知数个数。

---

## 3. 算例一：2 个方程 2 个未知数

### 3.1 方程组与矩阵

$$
\begin{cases}
2x - y = 0\\
-x + 2y = 3
\end{cases}
\quad\Longrightarrow\quad
A=\begin{bmatrix}2&-1\\-1&2\end{bmatrix},\quad
\mathbf{x}=\begin{bmatrix}x\\y\end{bmatrix},\quad
\mathbf{b}=\begin{bmatrix}0\\3\end{bmatrix}
$$

### 3.2 行图像

在 $xy$ 平面上，**每个方程画出满足它的全部点，构成一条直线**（"linear" 一词本来就含 line 的意思）。

- 方程 1：$2x-y=0$，即 $y=2x$。
  - 取 $y=0$ 得 $(0,0)$；取 $x=1$ 得 $y=2$，即点 $(1,2)$；取 $x=\tfrac12$ 得 $y=1$。
  - **含 $0$ 的右端项使直线过原点**。
- 方程 2：$-x+2y=3$。
  - 取 $y=0$ 得 $x=-3$，即点 $(-3,0)$；取 $x=-1$ 得 $y=1$，即点 $(-1,1)$。
  - 右端项不为 $0$，**这条直线不过原点**（否则 $0=3$，矛盾）。

两直线的交点就是解：

$$
y=2x \;\Rightarrow\; -x+2(2x)=3 \;\Rightarrow\; 3x=3 \;\Rightarrow\; x=1,\; y=2
$$

代回验算：$-1+4=3$，成立。故解为 $(x,y)=(1,2)$。

> **易错点**：$n=2$ 时行图像很好用；但「每个方程一条直线」只对二维成立。三维中每个方程是**一张平面**，不是直线。

### 3.3 列图像

把 $A\mathbf{x}=\mathbf{b}$ 按**列**展开：

$$
x\begin{bmatrix}2\\-1\end{bmatrix}+y\begin{bmatrix}-1\\2\end{bmatrix}=\begin{bmatrix}0\\3\end{bmatrix}
$$

即：寻找两列的**线性组合（linear combination）**，使其等于 $\mathbf{b}$。

取 $x=1,\;y=2$（上面已求出的解）：

$$
1\cdot\begin{bmatrix}2\\-1\end{bmatrix}+2\cdot\begin{bmatrix}-1\\2\end{bmatrix}
=\begin{bmatrix}2\\-1\end{bmatrix}+\begin{bmatrix}-2\\4\end{bmatrix}
=\begin{bmatrix}0\\3\end{bmatrix}=\mathbf{b}
$$

几何作图方式（**首尾相接**，向量加法）：

1. 从原点画出 $\mathbf{c}_1=\begin{bmatrix}2\\-1\end{bmatrix}$（向右 2、向下 1）；
2. 在 $\mathbf{c}_1$ 的终点接上 $\mathbf{c}_2=\begin{bmatrix}-1\\2\end{bmatrix}$（向左 1、向上 2）；
3. 由于 $y=2$，再接第二个 $\mathbf{c}_2$；
4. 终点落在 $\begin{bmatrix}0\\3\end{bmatrix}$，正好是 $\mathbf{b}$。

**更深的问题**：如果允许 $x,y$ 取**任意值**，所有组合能覆盖哪些 $\mathbf{b}$？

> 答案：**填满整个平面 $\mathbb{R}^2$**。因为两列不成比例（不共线），可以朝任意方向「张成」平面，于是**任何 $\mathbf{b}$ 都可达**。这就是后面所有「可解性」讨论的原型。

---

## 4. 算例二：3 个方程 3 个未知数

### 4.1 方程组与矩阵

$$
\begin{cases}
2x - y \phantom{+ 0z} = 0\\
-x + 2y - z = -1\\
\phantom{-x} -3y + 4z = 4
\end{cases}
\quad\Longrightarrow\quad
A=\begin{bmatrix}2&-1&0\\-1&2&-1\\0&-3&4\end{bmatrix},\quad
\mathbf{b}=\begin{bmatrix}0\\-1\\4\end{bmatrix}
$$

> **字幕校正**：原字幕此处的第三行转录为「-3z, -3ys … 4zs is, say, 4」，且第二行缺了加号。按各列系数反推可知第三行为 $-3y+4z=4$：第二列（$y$ 的系数）为 $(-1,2,-3)$，第三列（$z$ 的系数）为 $(0,-1,4)$，与 $\mathbf{b}=(0,-1,4)$ 完全自洽。

### 4.2 行图像

三维空间中 $x,y,z$ 三个坐标轴，**每一行给出三个未知数的一张平面**。以第二行为例，$-x+2y-z=-1$ 上的三个点：

| 取法 | 点 | 验算 |
|---|---|---|
| $y=z=0$ | $(1,0,0)$ | $-1+0-0=-1$ |
| $x=y=0$ | $(0,0,1)$ | $0+0-1=-1$ |
| $x=z=0$ | $(0,-\tfrac12,0)$ | $0-1-0=-1$ |

三点定一张平面。于是：

- 第一、二张平面相交于**一条直线**；
- 第三张平面与这条直线相交于**一个点**。

因为这三张平面「不平行、不特殊」，它们交于唯一一点，那个点就是解。

**解可以直接验证**：取 $(x,y,z)=(0,0,1)$，

$$
2\cdot0-0=0,\qquad -0+2\cdot0-1=-1,\qquad -3\cdot0+4\cdot1=4
$$

三个方程全部成立。

> **易错点（行图像的两个限制）**：
> 1. 三维的行图像已经很难画清楚，四维以上根本无法可视化；
> 2. 即便在三维，「三张平面」也可能因为平行或退化成交线而**不交于一点**——此时无唯一解。这也是为什么不能只依赖行图像来判断可解性。

### 4.3 列图像

$$
x\begin{bmatrix}2\\-1\\0\end{bmatrix}
+y\begin{bmatrix}-1\\2\\-3\end{bmatrix}
+z\begin{bmatrix}0\\-1\\4\end{bmatrix}
=\begin{bmatrix}0\\-1\\4\end{bmatrix}
$$

**注意这个特意构造的例子**：第三列 $\begin{bmatrix}0\\-1\\4\end{bmatrix}$ 恰好**就是** $\mathbf{b}$。因此一眼可见解为

$$
(x,y,z)=(0,0,1)
$$

即「取第三列一份、前两列都不取」。三个向量分别是 $\mathbb{R}^3$ 中的方向，问的是「用哪组系数组合出 $\mathbf{b}$」。

> **关键认知**：列图像里的解 $\mathbf{x}$ 是**组合系数**，不是坐标点；行图像里的解才是交点坐标。同一个 $(0,0,1)$，在两幅图里扮演的是不同角色。

### 4.4 换一个右端项

保持同一个 $A$，只改右端项。取 $\mathbf{b}=$ 第一列 $+$ 第二列：

$$
\mathbf{b}_{\text{new}}=\begin{bmatrix}2\\-1\\0\end{bmatrix}+\begin{bmatrix}-1\\2\\-3\end{bmatrix}=\begin{bmatrix}1\\1\\-3\end{bmatrix}
$$

**列图像**：同一个 $A$，同样的三列，只是换一个目标。解法立刻可见：

$$
1\cdot\mathbf{c}_1+1\cdot\mathbf{c}_2+0\cdot\mathbf{c}_3=\mathbf{b}_{\text{new}}
\;\Longrightarrow\;
(x,y,z)=(1,1,0)
$$

**行图像**：**三张新的平面**（因为右端项变了），它们现在交于 $(1,1,0)$。可见「同一个 $A$、不同的 $\mathbf{b}$」在两幅图里的表现方式完全不同。

---

## 5. 三种视角对照

| 视角 | 一次看什么 | 一行/一列的含义 | 解体现在哪里 | 可扩展性 |
|---|---|---|---|---|
| 行图像（row picture） | $A$ 的一行 | 一个方程 = $\mathbb{R}^n$ 中的一张超平面 | 所有超平面的公共交点（唯一解） | 仅 $n=2$ 直观，$n\ge3$ 迅速失效 |
| 列图像（column picture） | $A$ 的一列 | 一个未知数 = 该列向量的系数 | 凑出 $\mathbf{b}$ 的那组系数 $\mathbf{x}$ | 任意维数都可用，理论核心 |
| 矩阵形式（matrix form） | 整体 $A\mathbf{x}=\mathbf{b}$ | 紧凑记号，行与列同时存在 | 解为 $\mathbf{x}$（$A$ 可逆时唯一） | 与维数无关的抽象表述 |

---

## 6. 核心问题：对每个右端项都能解吗

这是本讲留下的最大问题，也是**代数问题**：

> 对任意右端项 $\mathbf{b}$，$A\mathbf{x}=\mathbf{b}$ 是否都有解？

**翻译成列图像的语言**（同一问题的另一种说法）：

> $A$ 各列的**全部线性组合**，能否填满整个三维空间 $\mathbb{R}^3$？

因为「求解 $A\mathbf{x}=\mathbf{b}$」与「用 $A$ 的列组合出 $\mathbf{b}$」是同一件事。

### 6.1 非奇异（可逆）情形

对上一节的 $A$，答案是**能**。这样的矩阵称为

- **非奇异矩阵（non-singular matrix）**
- **可逆矩阵（invertible matrix）**

这是最理想的一类矩阵：三列提供三个「独立的方向」，其组合铺满整个空间。

### 6.2 奇异情形：什么时候会失败

**失败的条件**：三列向量落在**同一张平面**内。

> 成因：三个向量共面时，它们的任何线性组合仍然落在这张平面里，**不可能跳出平面**。第三列提供不了任何新方向，全部组合只能覆盖那张平面。

具体例子：若**第三列恰好等于第一列加第二列**，则第三列完全冗余，得不到任何新东西。

后果：

- 只有**位于该平面内**的 $\mathbf{b}$ 有解；
- 平面之外的 $\mathbf{b}$ 永远取不到，**无解**；
- 这样的矩阵称为**奇异矩阵（singular matrix）**，**不可逆（not invertible）**。

> **补充（教材结论）**：当 $\mathbf{b}$ 落在这张平面内时方程确实有解，但由于列向量线性相关，解**不唯一**——例如 $\mathbf{c}_3=\mathbf{c}_1+\mathbf{c}_2$ 时，$\mathbf{c}_1+\mathbf{c}_2-\mathbf{c}_3=\mathbf{0}$ 说明 $\mathbf{x}=(1,1,-1)$ 是齐次方程 $A\mathbf{x}=\mathbf{0}$ 的非零解，于是解集是一条直线（无穷多解）。因此退化情形只有两种结局：**无解**或**无穷多解**，绝不会有唯一解。

对应的行图像语言：三张平面不再交于一点（例如其中两张平行，或交线落在第三张平面内）。

> **易错点**：「$n$ 个方程 $n$ 个未知数」**不代表**一定有唯一解。本讲的 $3\times3$ 例子同样是方程数等于未知数个数，它却可能在列共面时退化。**判断可解性要看列的组合能覆盖多大空间，而不是数方程的个数。**

### 6.3 推广到 9 维

若方程组有 9 个方程、9 个未知数：

- 有 9 个列向量，每个都是 $\mathbb{R}^9$ 中的向量；
- 问题变成：$\mathbb{R}^9$ 中 9 个向量的全部线性组合，能否覆盖整个 $\mathbb{R}^9$？
- 若矩阵是「随机」取的（例如用 MATLAB 的 `rand` 生成一个 $9\times9$ 矩阵），几乎必然是好的：非奇异、可逆。
- 但若第 9 列与第 8 列**相同**（或线性相关），它就不带来新方向，组合只能铺满 $\mathbb{R}^9$ 中的一个 **8 维超平面**，于是有些 $\mathbf{b}$ 取不到。

> **这就是线性代数的中心思想**：即使无法真正可视化高维空间，也应习惯「$n$ 个 $n$ 维向量的线性组合能覆盖多大空间」这种提问方式。

---

## 7. 矩阵乘向量：结果就是各列的线性组合

本讲最后回到矩阵形式，并给出**矩阵乘向量**的两种算法。以

$$
A=\begin{bmatrix}2&5\\1&3\end{bmatrix},\qquad \mathbf{x}=\begin{bmatrix}1\\2\end{bmatrix}
$$

为例。

### 7.1 列方法（推荐，也是本课程主视角）

$$
A\mathbf{x}=1\cdot\begin{bmatrix}2\\1\end{bmatrix}+2\cdot\begin{bmatrix}5\\3\end{bmatrix}
=\begin{bmatrix}2\\1\end{bmatrix}+\begin{bmatrix}10\\6\end{bmatrix}
=\begin{bmatrix}12\\7\end{bmatrix}
$$

即：按 $\mathbf{x}$ 的分量去**加权各列**再相加。

### 7.2 行方法（点积 / 数量积）

$$
A\mathbf{x}=\begin{bmatrix}2\cdot1+5\cdot2\\1\cdot1+3\cdot2\end{bmatrix}=\begin{bmatrix}12\\7\end{bmatrix}
$$

每一行与 $\mathbf{x}$ 作**点积（dot product）**。两种算法结果一致，不是巧合——它们互为定义的不同展开。

### 7.3 两种算法对照

| 算法 | 操作 | 得到的分量 | 何时用 |
|---|---|---|---|
| 列方法 | $\sum_j x_j\,\mathbf{c}_j$：$A$ 的列的线性组合 | 整体向量一次得到 | 理论推导、列空间、无解判定 |
| 行方法 | 每一行与 $\mathbf{x}$ 做点积 | 逐个分量得到 | 手算、编程实现、与方程逐行对应 |

### 7.4 结论

$$
\boxed{\;A\mathbf{x} \text{ 是 } A \text{ 的列的线性组合}\;}
$$

这是本讲最重要的一句话：**看到 $A\mathbf{x}=\mathbf{b}$，就想到「用 $A$ 的列凑出 $\mathbf{b}$」。**

顺带地，若 $\mathbf{b}=\begin{bmatrix}12\\7\end{bmatrix}$，则解自然是 $\mathbf{x}=\begin{bmatrix}1\\2\end{bmatrix}$——这正是上面算例的逆读。

---

## 8. 术语对照（本讲出现的专业名词）

| 中文（通行译名） | English | 说明 |
|---|---|---|
| 线性方程组 | system of linear equations | 本讲的根本问题 |
| 未知数 | unknowns | 个数决定 $A$ 的列数 |
| 系数矩阵 | coefficient matrix | $A$，按行排列系数 |
| 右端项 / 右端向量 | right-hand side | $\mathbf{b}$ |
| 向量 | vector | 有分量的一列数 |
| 行图像 | row picture | 一次看一个方程，解是交点/交线 |
| 列图像 | column picture | 一次看一列，解是组合系数 |
| 线性组合 | linear combination | 用系数加权若干向量再相加 |
| 矩阵形式 | matrix form | $A\mathbf{x}=\mathbf{b}$ |
| 平面 | plane | 三维中一个线性方程的解集 |
| 非奇异（的） | non-singular | 各列提供独立方向，$A\mathbf{x}=\mathbf{b}$ 对每个 $\mathbf{b}$ 可解 |
| 可逆（的） | invertible | 与「非奇异」同义 |
| 奇异（的） | singular | 列向量共面/相关，存在取不到的 $\mathbf{b}$ |
| 点积（数量积） | dot product | 行方法与 $\mathbf{x}$ 相乘的方式 |
| 消元法 | elimination | 下一讲的系统化求解方法 |
| 列空间 | column space | 后续讲次正式命名「所有线性组合构成的集合」 |
| 增广矩阵 | augmented matrix | 下一讲消元法中会用到 |

---

## 总结

1. **根本问题**：求解 $n$ 个方程、$n$ 个未知数的线性方程组，统一写成 $A\mathbf{x}=\mathbf{b}$。
2. **行图像**：每个方程在 $\mathbb{R}^n$ 中给出一个超平面（二维是直线，三维是平面），解是它们的公共交点；二维直观，三维以上迅速失效。
3. **列图像（最重要）**：求解 $A\mathbf{x}=\mathbf{b}$ 等价于用 $A$ 的各列线性组合出 $\mathbf{b}$，解 $\mathbf{x}$ 就是那组组合系数；这一视角不受维数限制。
4. **$2\times2$ 算例**：$A=\begin{bmatrix}2&-1\\-1&2\end{bmatrix}$、$\mathbf{b}=\begin{bmatrix}0\\3\end{bmatrix}$，两列组合得 $\mathbf{b}$ 的唯一系数是 $(x,y)=(1,2)$。
5. **$3\times3$ 算例**：$A=\begin{bmatrix}2&-1&0\\-1&2&-1\\0&-3&4\end{bmatrix}$、$\mathbf{b}=\begin{bmatrix}0\\-1\\4\end{bmatrix}$，因第三列恰等于 $\mathbf{b}$，解为 $(0,0,1)$；换 $\mathbf{b}=\mathbf{c}_1+\mathbf{c}_2=\begin{bmatrix}1\\1\\-3\end{bmatrix}$ 则解为 $(1,1,0)$。
6. **可解性判据**：$A$ 的列能否铺满整个空间。列的线性组合铺满全空间 ⇔ 非奇异（可逆），对每个 $\mathbf{b}$ 可解；若各列共面（例如 $\mathbf{c}_3=\mathbf{c}_1+\mathbf{c}_2$），则平面外的 $\mathbf{b}$ **无解**、平面内的 $\mathbf{b}$ 有**无穷多解**——绝不会唯一。方程数等于未知数个数**不保证**唯一解。
7. **核心公式**：$A\mathbf{x}=\mathbf{b}$ 中 $A\mathbf{x}$ 是 $A$ 的列的线性组合；下一讲用**消元法**系统化地求解并判断何时无解。

## 参考资料

- [MIT OCW 18.06 Linear Algebra 课程主页（视频、讲义、习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/)
- [本讲页面：Lecture 1: The Geometry of Linear Equations（视频、讲义与习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-1-the-geometry-of-linear-equations/)
- [MIT 18.06 官方课程网站（web.mit.edu/18.06，本讲中提到的课程主页）](https://web.mit.edu/18.06/www/)
- [MIT 18.06 课程要点摘要（mitmath/1806 仓库）](https://github.com/mitmath/1806)
- [Gilbert Strang《Introduction to Linear Algebra》教材官网（目录与配套材料）](https://math.mit.edu/~gs/linearalgebra/)
- [MIT 18.06 线性代数中文笔记（ApacheCN 翻译整理，供中文术语对照）](https://github.com/apachecn/mit-18.06-linalg-notes)

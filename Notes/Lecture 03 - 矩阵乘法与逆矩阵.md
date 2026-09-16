---
title: "MIT 18.06 Lecture 03: Multiplication and Inverse Matrices"
course: "MIT 18.06 Linear Algebra (Gilbert Strang)"
book: "Introduction to Linear Algebra, Gilbert Strang"
tags: [线性代数, 18.06, MIT, 矩阵乘法, 逆矩阵, 高斯-若尔当消元, 奇异矩阵, 分块矩阵, 初等变换]
---

# Lecture 03 — 矩阵乘法与逆矩阵

> **本讲两个主题**
> 1. **矩阵乘法**：有多种算法（逐元素、整列、整行、列乘行、分块），结果完全相同，各有适用场合。
> 2. **逆矩阵（inverse matrix）**：方阵何时可逆（invertible）？如何求？——**高斯-若尔当消元法（Gauss-Jordan elimination）**。

---

## 1. 尺寸规则与定义

设 $A$ 是 $m\times n$ 矩阵，$B$ 是 $n\times p$ 矩阵，则乘积 $C=AB$ 是 $m\times p$ 矩阵。

> **唯一硬性约束**：左边矩阵的**列数**必须等于右边矩阵的**行数**（$n$ 对 $n$）。方阵相乘则必须同阶。

$C$ 的第 $i$ 行第 $j$ 列元素（逐元素定义）：

$$c_{ij}=\sum_{k=1}^{n} a_{ik}b_{kj}$$

下标顺序永远是「先行号、后列号」：$c_{34}$ 取自 $A$ 的第 3 行与 $B$ 的第 4 列的内积：

$$c_{34}=a_{31}b_{14}+a_{32}b_{24}+\cdots+a_{3n}b_{n4}$$

> 这是「标准算法」，大家（包括 Strang 本人）平时都这么算。但**它不是唯一视角**，另外几种视角在后续课程（列空间、行空间、消元、分块）中反复用到。

---

## 2. 矩阵乘法的五种视角

| 视角 | 公式 | 记忆要点 | 典型用途 |
|---|---|---|---|
| ① 逐元素内积（行 × 列） | $c_{ij}=\sum_k a_{ik}b_{kj}$ | 行点乘列，得到**一个数** | 手算、编程实现 |
| ② 整列（矩阵 × 向量） | $C_{:,j}=A\,B_{:,j}$ | $C$ 的每一列是 **$A$ 的列的组合**，系数取自 $B$ 的该列 | 列空间、$Ax=b$ |
| ③ 整行（向量 × 矩阵） | $C_{i,:}=A_{i,:}\,B$ | $C$ 的每一行是 **$B$ 的行的组合**，系数取自 $A$ 的该行 | 行空间、消元 |
| ④ 列乘行（外积） | $AB=\sum_{k=1}^{n}(\text{col}_k A)(\text{row}_k B)$ | 列（$m\times1$）× 行（$1\times p$）= **一个满尺寸 $m\times p$ 矩阵** | 低秩分解、$A=uv^{\mathsf T}$ |
| ⑤ 分块乘法 | 见 §2.5 | **把块当元素**照样相乘相加 | 大规模计算、证明 |

### 2.1 视角① 逐元素内积

只对单个元素而言，$c_{34}$ 是 $A$ 第 3 行与 $B$ 第 4 列的内积。它的意义在于「可计算」；缺点是完全看不出整体结构。

### 2.2 视角② 按整列看

把 $B$ 看成 $p$ 个列向量并排摆放，则

$$A\begin{bmatrix} \mid & \mid & & \mid \\ b_1 & b_2 & \cdots & b_p \\ \mid & \mid & & \mid \end{bmatrix}
=\begin{bmatrix} \mid & \mid & & \mid \\ Ab_1 & Ab_2 & \cdots & Ab_p \\ \mid & \mid & & \mid \end{bmatrix}$$

**为什么成立**：$Ab_j$ 的定义就是 $A$ 各列按 $b_j$ 的分量做线性组合；而视角①中 $c_{ij}$ 只用到了 $A$ 的第 $i$ 行和 $b_j$，与 $B$ 的其他列无关。

> **结论**：$AB$ 的每一列都是 $A$ 的列向量的线性组合，系数就是 $B$ 对应列的元素。$A$ 的列向量有 $m$ 个分量，$C$ 的列向量也有 $m$ 个分量——**列数可以变，行数不变**。

**算例**：取

$$A=\begin{bmatrix}2&7\\3&8\\4&9\end{bmatrix}\ (3\times2),\qquad
B=\begin{bmatrix}1&6\\0&0\end{bmatrix}\ (2\times2)$$

$B$ 的第 1 列是 $\begin{bmatrix}1\\0\end{bmatrix}$，第 2 列是 $\begin{bmatrix}6\\0\end{bmatrix}$。

- $C$ 第 1 列 $=1\cdot\begin{bmatrix}2\\3\\4\end{bmatrix}+0\cdot\begin{bmatrix}7\\8\\9\end{bmatrix}=\begin{bmatrix}2\\3\\4\end{bmatrix}$
- $C$ 第 2 列 $=6\cdot\begin{bmatrix}2\\3\\4\end{bmatrix}+0\cdot\begin{bmatrix}7\\8\\9\end{bmatrix}=\begin{bmatrix}12\\18\\24\end{bmatrix}$

$$AB=\begin{bmatrix}2&12\\3&18\\4&24\end{bmatrix}$$

### 2.3 视角③ 按整行看

$$\begin{bmatrix}- & a_1 & -\\- & a_2 & -\\ & \vdots & \end{bmatrix}B
=\begin{bmatrix}- & a_1B & -\\- & a_2B & -\\ & \vdots & \end{bmatrix}$$

**结论**：$AB$ 的每一行都是 $B$ 的行的线性组合，系数取自 $A$ 的对应行。左乘 $A$ 的作用是「把 $B$ 的行混合起来」。

**算例**（同一组 $A,B$）：$A$ 的第 1 行是 $\begin{bmatrix}2&7\end{bmatrix}$，$B$ 的两行是 $\begin{bmatrix}1&6\end{bmatrix}$ 与 $\begin{bmatrix}0&0\end{bmatrix}$，故

$$C \text{ 第 1 行}=2\begin{bmatrix}1&6\end{bmatrix}+7\begin{bmatrix}0&0\end{bmatrix}=\begin{bmatrix}2&12\end{bmatrix}$$

与 §2.2 结果一致。

> **对称记忆法**：**右乘按列组合 $A$ 的列，左乘按行组合 $B$ 的行。** 这是本讲最值得背下来的一句。

### 2.4 视角④ 列乘行（外积之和）

**行 × 列 → 一个数；列 × 行 → 一个完整矩阵。**

$$\begin{bmatrix}2\\3\\4\end{bmatrix}\begin{bmatrix}1&6\end{bmatrix}
=\begin{bmatrix}2&12\\3&18\\4&24\end{bmatrix}$$

这个矩阵极其特殊：三行都是 $\begin{bmatrix}1&6\end{bmatrix}$ 的倍数，两列都是 $\begin{bmatrix}2\\3\\4\end{bmatrix}$ 的倍数。用后续课程的语言说：

- **行空间（row space）** 就是过向量 $\begin{bmatrix}1&6\end{bmatrix}$ 的那条直线；
- **列空间（column space）** 就是过向量 $\begin{bmatrix}2\\3\\4\end{bmatrix}$ 的那条直线。

因此它是一个「最小」的矩阵——秩为 1。这也解释了 §2.2、§2.3 的结论：列组合只有一个向量可选时，就退化为「倍数」。

**一般公式**：

$$AB=\sum_{k=1}^{n}(\text{第 }k\text{ 列 of }A)(\text{第 }k\text{ 行 of }B)$$

**跟算**：用 §2.2 的同一组矩阵。

$$AB=\begin{bmatrix}2\\3\\4\end{bmatrix}\begin{bmatrix}1&6\end{bmatrix}+\begin{bmatrix}7\\8\\9\end{bmatrix}\begin{bmatrix}0&0\end{bmatrix}
=\begin{bmatrix}2&12\\3&18\\4&24\end{bmatrix}+\begin{bmatrix}0&0\\0&0\\0&0\end{bmatrix}
=\begin{bmatrix}2&12\\3&18\\4&24\end{bmatrix}$$

结果与视角②③完全一致。

### 2.5 视角⑤ 分块乘法（block multiplication）

把 $A$、$B$ 各切成四块（例如 $20\times20$ 的矩阵切成 $10\times10$ 的块）：

$$A=\begin{bmatrix}A_1&A_2\\A_3&A_4\end{bmatrix},\qquad B=\begin{bmatrix}B_1&B_2\\B_3&B_4\end{bmatrix}$$

则

$$AB=\begin{bmatrix}A_1B_1+A_2B_3 & A_1B_2+A_2B_4\\ A_3B_1+A_4B_3 & A_3B_2+A_4B_4\end{bmatrix}$$

**规则**：把每个块当成一个「元素」，完全照搬 $2\times2$ 矩阵的乘法公式。左上块的第一个块乘积就是 $A_1B_1$，再加上 $A_2B_3$。

> **为什么能这样**：逐项展开后，五种算法做的是**同一批乘法**，只是分组方式不同。分块唯一的约束是**分块的尺寸必须相互匹配**（相乘的块之间列数 = 行数），块本身不必同形。

---

## 3. 逆矩阵

### 3.1 定义与基本事实

设 $A$ 为**方阵**。若存在矩阵（记作 $A^{-1}$）使

$$A^{-1}A=I \quad\text{且}\quad AA^{-1}=I$$

则称 $A$ **可逆（invertible）**，也称**非奇异（nonsingular）**或**满秩（full rank）**；$A^{-1}$ 称为 $A$ 的**逆矩阵（inverse matrix）**。不存在这样的矩阵时，称 $A$ **不可逆**或**奇异矩阵（singular matrix）**。

> **对每个方阵，最重要的问题就是：它可逆吗？** 这是本讲的核心问题，也是整门课反复回到的问题。

**关于左逆与右逆**（Strang 特别强调）：

- 对于**方阵**：只要存在矩阵 $C$ 使 $CA=I$，那么同一个 $C$ 从右边乘也得到 $I$（即 $AC=I$），所以**左逆就是右逆**。（这一点并不容易证明，但成立。）
- 对于**长方形矩阵**：可能出现「左逆」与「右逆」不是同一个矩阵的情况——事实上形状本身就不允许同一个矩阵两边都成立。

### 3.2 奇异矩阵：一个 $2\times2$ 反例

$$A=\begin{bmatrix}1&3\\2&6\end{bmatrix}$$

**判据一（行列式，本课尚未讲）**：$\det A = 1\cdot6-3\cdot2=0$。

**判据二（列的角度）**：两列分别是 $\begin{bmatrix}1\\2\end{bmatrix}$ 与 $\begin{bmatrix}3\\6\end{bmatrix}$，且

$$\begin{bmatrix}3\\6\end{bmatrix}=3\begin{bmatrix}1\\2\end{bmatrix}$$

两列在同一条直线上。由视角②，$AB$ 的列只能是 $A$ 的列的线性组合，永远落在这条直线上；而 $I$ 的列 $\begin{bmatrix}1\\0\end{bmatrix}$ 不在这条直线上。**所以 $A$ 不可能有逆。**

**判据三（零空间角度，最本质）**：存在**非零**向量 $x$ 使 $Ax=0$。此处

$$x=\begin{bmatrix}3\\-1\end{bmatrix},\qquad
Ax=\begin{bmatrix}1&3\\2&6\end{bmatrix}\begin{bmatrix}3\\-1\end{bmatrix}
=\begin{bmatrix}3-3\\6-6\end{bmatrix}=\begin{bmatrix}0\\0\end{bmatrix}$$

（$x=\begin{bmatrix}0\\0\end{bmatrix}$ 不算，它永远成立，没有信息量。）

> **灾难性推理**：若 $A^{-1}$ 存在，对 $Ax=0$ 两边左乘 $A^{-1}$，得
> $$x=A^{-1}0=0$$
> 与 $x\neq 0$ 矛盾。可见 $A^{-1}$ 根本不可能存在——**零永远被零带到零，逆矩阵无法把 0 还原成 $x$**。

**结论**：非奇异（可逆）$\iff$ 不存在非零 $x$ 使 $Ax=0$ $\iff$ $A$ 的各列线性无关。奇异矩阵必存在列的组合等于零列。

### 3.3 可逆矩阵与逆矩阵的求法（分为两个方程组）

取一个确实可逆的例子：

$$A=\begin{bmatrix}1&3\\2&7\end{bmatrix},\qquad \det A=1\cdot7-3\cdot2=1\neq0$$

两个列向量指向不同方向，因此可以组合出任意向量。

记 $A^{-1}=\begin{bmatrix}a&c\\b&d\end{bmatrix}$。由 $AA^{-1}=I$ 按**列**拆开（视角②）：

$$A\begin{bmatrix}a\\b\end{bmatrix}=\begin{bmatrix}1\\0\end{bmatrix},\qquad
A\begin{bmatrix}c\\d\end{bmatrix}=\begin{bmatrix}0\\1\end{bmatrix}$$

即 $A\times(A^{-1}\text{的第 }j\text{ 列})=I$ 的第 $j$ 列。**求逆矩阵 = 用同一个系数矩阵 $A$、同时解 $n$ 个右端项不同的方程组。**

### 3.4 高斯-若尔当消元：$[A\mid I]\to[I\mid A^{-1}]$

做法：把 $I$ 作为额外的列拼到 $A$ 右边得到**增广矩阵（augmented matrix）**，对左半部分做**初等行变换（elementary row operation）**化为 $I$，则右半部分自动变成 $A^{-1}$。

$$\begin{bmatrix}A\mid I\end{bmatrix}=\begin{bmatrix}1&3&1&0\\2&7&0&1\end{bmatrix}$$

**第 1 步（向下消元，把 $a_{21}$ 消成 0）**：第 2 行减去第 1 行的 2 倍

$$\begin{bmatrix}1&3&1&0\\0&1&-2&1\end{bmatrix}$$

到此为止左半部分是**上三角**——普通高斯消元（Gauss elimination）在这里就会停下。

> **Jordan 的贡献**：不要停，继续**向上消元**。

**第 2 步（向上消元，把 $a_{12}=3$ 消成 0）**：第 1 行减去第 2 行的 3 倍

$$\begin{bmatrix}1&0&7&-3\\0&1&-2&1\end{bmatrix}$$

于是

$$A^{-1}=\begin{bmatrix}7&-3\\-2&1\end{bmatrix}$$

**验算**（务必自己乘一遍）：

$$A^{-1}A=\begin{bmatrix}7&-3\\-2&1\end{bmatrix}\begin{bmatrix}1&3\\2&7\end{bmatrix}
=\begin{bmatrix}7-6&21-21\\-2+2&-6+7\end{bmatrix}
=\begin{bmatrix}1&0\\0&1\end{bmatrix}=I$$

### 3.5 为什么 $[A\mid I]$ 的右半部分会变成 $A^{-1}$

每一步消元都等价于左乘一个**初等矩阵（消元矩阵）$E$**。本例中有两个：

$$E_1=\begin{bmatrix}1&0\\-2&1\end{bmatrix}\ (\text{第 2 行}-2\times\text{第 1 行}),\qquad
E_2=\begin{bmatrix}1&-3\\0&1\end{bmatrix}\ (\text{第 1 行}-3\times\text{第 2 行})$$

整个过程相当于左乘 $E=E_2E_1$，作用在增广矩阵上：

$$E[A\mid I]=[EA\mid E]$$

而消元的目标是让左半部分变成 $I$，即 $EA=I$。由 $EA=I$ 立刻可知

$$E=A^{-1}$$

所以右半部分 $EI=E=A^{-1}$。

> **一句话**：$[A\mid I]\to[I\mid A^{-1}]$ 成立的唯一理由是——**把 $A$ 变成 $I$ 的那个变换，按定义就是 $A^{-1}$**。

**用矩阵乘法核对**（本讲两次消元的 $E$ 之积确实等于上一步求出的 $A^{-1}$）：

$$E=E_2E_1=\begin{bmatrix}1&-3\\0&1\end{bmatrix}\begin{bmatrix}1&0\\-2&1\end{bmatrix}
=\begin{bmatrix}1+6&-3\\-2&1\end{bmatrix}
=\begin{bmatrix}7&-3\\-2&1\end{bmatrix}=A^{-1}$$

---

## 4. 可逆性判据对照

| 判据 | 可逆（非奇异） | 不可逆（奇异） |
|---|---|---|
| 行列式 | $\det A\neq0$ | $\det A=0$ |
| 列向量 | 各列**线性无关**，指向不同方向 | 存在列的组合 = 零列（列共线/退化） |
| 零空间 | 只有 $x=0$ 满足 $Ax=0$ | 存在**非零** $x$ 使 $Ax=0$ |
| 能否解出 $I$ | 由 $A$ 的列能组合出 $I$ 的每一列 | $I$ 的列不在 $A$ 的列空间中，永远组合不出 |
| 消元结果 | $[A\mid I]\to[I\mid A^{-1}]$，主元全非零 | 消元出现**零主元**，左半部分无法化为 $I$ |

> **两点提醒**
> 1. 求逆矩阵时**必须用同一个 $A$ 解 $n$ 个方程组**，右端项分别是 $I$ 的各列——这正是「$n$ 个右端项」的含义，也是 Gauss-Jordan 相比「逐个解」高效的原因。
> 2. 对长方形矩阵，左逆与右逆不是同一个矩阵（形状也不允许），**不要把方阵的结论想当然地推广**。

---

## 5. 补充：逆矩阵的乘积法则（教材内容，本讲未展开）

若 $A$、$B$ 都是 $n$ 阶可逆矩阵，则 $AB$ 也可逆，且

$$(AB)^{-1}=B^{-1}A^{-1}$$

**为什么顺序反了**：要「撤销」先 $B$ 后 $A$ 的两次连续变换，必须**先撤销后做的那个**——就像穿脱鞋袜：穿时先穿袜再穿鞋，脱时先脱鞋再脱袜，顺序必然相反。验证只需直接相乘：

$$(AB)(B^{-1}A^{-1})=A(BB^{-1})A^{-1}=AIA^{-1}=AA^{-1}=I$$

> 一般地：$(A_1A_2\cdots A_k)^{-1}=A_k^{-1}\cdots A_2^{-1}A_1^{-1}$。**前提是所有 $A_i$ 都可逆且为同阶方阵**；乘积中只要有一个奇异，整体就不可逆。

---

## 总结

1. **矩阵乘法只有一个定义，但有五种等价视角**：逐元素内积、按整列、按整行、列乘行之和、分块。结果相同，用途不同——列视角通向列空间与 $Ax=b$，行视角通向行空间与消元，列乘行通向秩 1 矩阵与低秩分解，分块视角是大规模计算的基础。
2. **两句话抓住列与行**：右乘 $B$，把 $A$ 的列按 $B$ 的列系数组合；左乘 $A$，把 $B$ 的行按 $A$ 的行系数组合。
3. **行 × 列得一个数，列 × 行得一个满尺寸矩阵**。后者（外积）的所有列共线、所有行共线，秩为 1，是理解秩与行/列空间的最小例子。
4. **方阵最重要的问题是「可逆吗」**。可逆判据等价于：$\det A\neq0$、各列线性无关、$Ax=0$ 只有零解、消元不出现零主元。
5. **奇异的本质是「有非零 $x$ 被压成 0」**，因此 $A^{-1}$ 无从还原，逆矩阵必然不存在。
6. **求逆矩阵 = 解 $n$ 个同系数、不同右端项的方程组**；高斯-若尔当消元把它压缩成对增广矩阵 $[A\mid I]$ 的统一行变换，化为 $[I\mid A^{-1}]$。其正确性来自 $EA=I\Rightarrow E=A^{-1}$。
7. **下一讲预告**：把这些消元矩阵 $E$ 组织起来，就得到 $A=LU$ 分解（Lecture 04）。

---

## 参考资料

- [MIT OCW 18.06 Linear Algebra 课程主页（视频、讲义、习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/)
- [本讲页面：Lecture 3: Multiplication and Inverse Matrices（视频、讲义与习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-3-multiplication-and-inverse-matrices/)
- [下一讲页面：Lecture 4: Factorization into A = LU（消元矩阵的乘积即 LU 分解）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-4-factorization-into-a-lu/)
- [MIT 18.06 官方课程网站（web.mit.edu/18.06）：按主题索引的视频与讲义](https://web.mit.edu/18.06/www/)
- [MIT 18.06 课程要点摘要（mitmath/1806 仓库）](https://github.com/mitmath/1806)
- [Gilbert Strang《Introduction to Linear Algebra》教材官网（目录与配套材料）](https://math.mit.edu/~gs/linearalgebra/)
- [MIT 18.06 线性代数中文笔记（ApacheCN 翻译整理，供中文术语对照）](https://github.com/apachecn/mit-18.06-linalg-notes)

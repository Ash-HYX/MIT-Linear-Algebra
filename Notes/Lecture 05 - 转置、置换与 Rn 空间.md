---
title: "MIT 18.06 Lecture 05: Transposes, Permutations, Spaces R^n"
course: "MIT 18.06 Linear Algebra (Gilbert Strang)"
book: "Introduction to Linear Algebra, Gilbert Strang"
tags: [线性代数, 18.06, MIT, 转置, 向量空间, 置换矩阵, 对称矩阵, 子空间, 列空间, 消元法]
---

# Lecture 05 — 转置、置换与 R^n 空间

> 本讲分三段：先用**置换矩阵（permutation matrix）**补完消元法最后一块拼图，得到最一般的消元公式 $PA=LU$；再讲**转置（transpose）**与**对称矩阵（symmetric matrix）**；最后进入全书真正的主线——**向量空间（vector space）**与**子空间（subspace）**。

---

## 1. 置换矩阵与 PA = LU

### 1.1 为什么需要行交换

消元法分解 $A=LU$ 有一个隐含前提：**主元（pivot）位置上不出现零**。

- 如果主元位置出现 $0$，就把它下面的某一行换上来，得到一个非零主元，再继续消元。
- 换一次不够，可能要换两次甚至更多次。
- 只要允许「换行」，任意可逆矩阵都能走完消元流程。

因此必须把行交换这件事写进公式里。

### 1.2 一般情形：PA = LU

$$A = LU \quad\Longrightarrow\quad PA = LU$$

| 量 | 含义 |
|---|---|
| $P$ | 置换矩阵：记录**行交换**，把行调整成「主元不出零」的顺序 |
| $L$ | 下三角矩阵，对角元全为 $1$，下方是乘数（multipliers） |
| $U$ | 上三角矩阵，主对角线下方全为 $0$ |

> **易错点**：$A=LU$ 是「不需要换行」的特殊情形，此时 $P=I$。千万不要把 $A=LU$ 当成普遍公式——普遍公式是 $PA=LU$。绝大多数可逆矩阵不需要 $P$，但总有一小部分必须换行。

### 1.3 关于置换矩阵的三条硬事实

1. **$P$ 是单位矩阵重排行之后得到的矩阵**：每一行、每一列恰好有一个 $1$，其余为 $0$。「不换行」也算一种置换，即 $P=I$。
2. **$n\times n$ 置换矩阵共有 $n!$ 个**：
   $$n(n-1)(n-2)\cdots 3\cdot 2\cdot 1 = n!$$
   例：$4\times4$ 有 $4!=24$ 个；$5\times5$ 有 $5!=120$ 个，靠逐个列举不现实。
3. **每个置换矩阵都可逆，且逆矩阵就是它的转置**：
   $$P^TP = I \quad\Longleftrightarrow\quad P^{-1} = P^{T}$$
   直觉：转置把「$1$ 的位置」沿主对角线镜像，正好把行还原回原来的顺序。

**算例（交换两行）**

$$P_{12}=\begin{bmatrix}0&1&0\\1&0&0\\0&0&1\end{bmatrix},\qquad
P_{12}^{T}=\begin{bmatrix}0&1&0\\1&0&0\\0&0&1\end{bmatrix}=P_{12},\qquad
P_{12}^{T}P_{12}=I$$

**算例（三轮换 $1\to2\to3\to1$）**

$$P=\begin{bmatrix}0&0&1\\1&0&0\\0&1&0\end{bmatrix},\qquad
P^{T}=\begin{bmatrix}0&1&0\\0&0&1\\1&0&0\end{bmatrix},\qquad
P^{T}P=\begin{bmatrix}1&0&0\\0&1&0\\0&0&1\end{bmatrix}=I$$

此时 $P$ 的逆**不等于** $P$ 自己，但一定等于 $P^{T}$。

> **后续**：满足 $M^{T}M=I$（即列向量标准正交）的矩阵族比置换矩阵大得多，称为**正交矩阵（orthogonal matrix）**；置换矩阵只是其中一个特殊子集。

### 1.4 数值视角：为什么 Matlab 会「多此一举」地换行

- 人判断主元只看「是不是零」；Matlab 还判断「**主元是不是足够大**」。
- 原因：**接近零的主元在数值上很危险**（会放大舍入误差，导致结果不可信）。
- 所以 Matlab 解方程组时会做一些「代数上没必要、但精度上必要」的行交换——这叫**选主元（pivoting）**。
- 而课程在此只讨论精确代数，所以结论是：**能不用 $P$ 就不用，但必须允许用。**

---

## 2. 转置与对称矩阵

### 2.1 转置的定义

$$(A^{T})_{ij} = A_{ji}$$

即**沿主对角线翻转**：行号与列号互换。$m\times n$ 矩阵的转置是 $n\times m$。

**算例（黑板原题）**

$$A=\begin{bmatrix}1&2&4\\3&3&1\end{bmatrix}_{2\times3}
\qquad\Longrightarrow\qquad
A^{T}=\begin{bmatrix}1&3\\2&3\\4&1\end{bmatrix}_{3\times2}$$

- $A$ 的第一行 $[1,2,4]$ 变成 $A^{T}$ 的第一列；
- $A$ 的每一列变成 $A^{T}$ 的对应行；
- 形状交换：矮宽的变成高瘦的。

> **字幕校勘**：自动转录中把 $A$ 口述成「three by two / 变宽」，与黑板上写出的矩阵形状矛盾。以写出来的矩阵为准：$A$ 是 $2\times3$，$A^{T}$ 是 $3\times2$。后面的 $R^{T}R$ 算例（得到 $3\times3$）也印证了这一点。

### 2.2 转置的运算律

转置乘积**顺序反转**，和求逆的规则完全一样：

$$(AB)^{T} = B^{T}A^{T}, \qquad (A^{T})^{T} = A$$

### 2.3 对称矩阵（symmetric matrix）

**定义**：转置后不变，即

$$S^{T} = S \quad\Longleftrightarrow\quad S_{ij} = S_{ji}$$

性质对照（可与「转置给出逆」的置换矩阵对照记忆）：

| 矩阵类型 | 满足的关系 | 直观 | 识别难度 |
|---|---|---|---|
| 置换矩阵 $P$ | $P^{T}P=I$，即 $P^{T}=P^{-1}$ | 转置给出**逆** | 需要做乘法才看得出 |
| 对称矩阵 $S$ | $S^{T}=S$ | 转置给出**自身** | 看一眼就认出：沿对角线镜像相同 |

**例（黑板形式）**：对角元任取，对角线两侧的数字必须成对镜像出现（黑板上用了 $1,7,9$ 这几个数）：

$$S=\begin{bmatrix}d_1&1&7\\1&d_2&9\\7&9&d_3\end{bmatrix},\qquad S^{T}=S$$

$S$ 中的 $1,7,9$ 同时出现在 $(1,2)/(2,1)$、$(1,3)/(3,1)$、$(2,3)/(3,2)$ 三个对称位置上——这就是对称的来源。

---

## 3. 核心结论：R^T R 一定是对称矩阵

**结论前置**：对任意矩阵 $R$（不必是方阵），$R^{T}R$ **永远是对称矩阵**，且是方阵。

**为什么**——证明只用了两行转置运算律：

$$(R^{T}R)^{T} = R^{T}(R^{T})^{T} = R^{T}R$$

转置之后没有变化，按定义它就是对称的。这也是实际应用中大量对称矩阵的来源。

**算例（跟算，黑板数字）**

取

$$R=\begin{bmatrix}1&2&4\\3&3&1\end{bmatrix},\qquad
R^{T}=\begin{bmatrix}1&3\\2&3\\4&1\end{bmatrix}$$

先算 $R^{T}R$ 的第一行（$R^{T}$ 的第一行 $[1,3]$ 分别点乘 $R$ 的三列）：

| 元素 | 计算 | 结果 |
|---|---|---|
| $(1,1)$ | $1\cdot1 + 3\cdot3$ | $10$ |
| $(1,2)$ | $1\cdot2 + 3\cdot3 = 2+9$ | $11$ |
| $(1,3)$ | $1\cdot4 + 3\cdot1 = 4+3$ | $7$ |

对角线与其余元素同理，最终得到：

$$R^{T}R=\begin{bmatrix}10&11&7\\11&13&11\\7&11&17\end{bmatrix}$$

- $(2,2)=2^{2}+3^{2}=13$，$(2,3)=2\cdot4+3\cdot1=11$，$(3,3)=4^{2}+1^{2}=17$；
- $11$ 在 $(1,2)$ 与 $(2,1)$ 出现两次，$7$ 在 $(1,3)$ 与 $(3,1)$ 出现两次——**这不是巧合，而是 $(AB)^{T}=B^{T}A^{T}$ 的必然结果**。

---

## 4. 向量空间与子空间

### 4.1 两种运算 = 线性组合

线性代数对向量只做两件事：

1. **加法（addition）**：$v+w$
2. **数乘（scalar multiplication）**：$cv$（$c$ 是实数，叫标量 scalar）

两者合起来就是**线性组合（linear combination）**。一个「向量的空间」必须允许对其中任意向量做这两件事，并且结果**不跑出去**。

> **术语**：这种「做完运算仍留在集合内」的性质叫**封闭（closed）**——加法封闭（closed under addition）、数乘封闭（closed under scalar multiplication）。

书本还列出了加法与数乘必须满足的八条规则，但这些规则在 $\mathbb{R}^{n}$ 中从不构成问题；**真正需要检验的永远是「封闭不封闭」**。

### 4.2 主要例子：R^n

$$\mathbb{R}^{n} = \{\text{所有含 } n \text{ 个实分量的向量}\}$$

- $\mathbb{R}^{2}$：所有二维实向量，几何上就是 $xy$ 平面。例：$\begin{bmatrix}3\\2\end{bmatrix}$、$\begin{bmatrix}0\\0\end{bmatrix}$、$\begin{bmatrix}\pi\\e\end{bmatrix}$。
- $\mathbb{R}^{3}$：所有三维实向量。
- $\mathbb{R}^{n}$：$n$ 个实分量。默认写成**列向量**。

**算例（逐分量相加）**

$$\begin{bmatrix}3\\2\end{bmatrix}+\begin{bmatrix}0\\0\end{bmatrix}=\begin{bmatrix}3\\2\end{bmatrix},\qquad
\begin{bmatrix}3\\2\end{bmatrix}+\begin{bmatrix}\pi\\e\end{bmatrix}=\begin{bmatrix}3+\pi\\2+e\end{bmatrix},\qquad
\begin{bmatrix}3\\2\end{bmatrix}+\begin{bmatrix}-3\\-2\end{bmatrix}=\begin{bmatrix}0\\0\end{bmatrix}$$

> **易错点**：$\begin{bmatrix}3\\2\\0\end{bmatrix}$ 属于 $\mathbb{R}^{3}$ 而不是 $\mathbb{R}^{2}$——分量个数说了算，某个分量恰好是 $0$ 并不降维。

### 4.3 反例：第一象限不是向量空间

取 $xy$ 平面中「两个分量都非负」的那四分之一：

$$Q=\left\{\begin{bmatrix}x\\y\end{bmatrix}\ \middle|\ x\ge0,\ y\ge0\right\}$$

- 加法没问题：$\begin{bmatrix}3\\2\end{bmatrix}+\begin{bmatrix}5\\6\end{bmatrix}=\begin{bmatrix}8\\8\end{bmatrix}$ 仍在 $Q$ 内；
- 数乘出问题：取标量 $-5$，
  $$-5\begin{bmatrix}3\\2\end{bmatrix}=\begin{bmatrix}-15\\-10\end{bmatrix}\notin Q$$
- 所以 $Q$ **不是**向量空间，原因是**对数乘不封闭**。

### 4.4 为什么每个向量空间都必须含零向量

零向量（zero vector）在空间里是**不可缺**的，理由有两条，都来自封闭性：

1. 允许用标量 $0$ 做数乘：$0\cdot v = 0$，结果必须留在空间内；
2. 允许把 $v$ 与它的相反向量 $-v$ 相加：$v+(-v)=0$，结果必须留在空间内。

> 因此，若某个集合在「取反/乘以 0」之后不含零向量，它就一定不是向量空间。这是检验子空间时**最快的一票否决**。

### 4.5 子空间（subspace）

**定义**：子空间是**某个向量空间的子集**，且它**自身仍是一个向量空间**——即对加法和数乘封闭。

关键推论：一旦某个非零向量 $v$ 属于子空间，则 $v$ 的**全部倍数 $cv$** 都必须在内，也就是**整条过原点的直线**必须在内。

**R^2 的全部子空间**

| 子空间 | 直观维数 | 说明 |
|---|---|---|
| 整个 $\mathbb{R}^{2}$ | 2 | 最大的子空间（空间本身算自己的子空间） |
| 过原点的任意直线 | 1 | **必须过原点**；不过原点的直线不是子空间（乘以 $0$ 就跑出去了） |
| 只含零向量的集合 $\{0\}$ | 0 | 最小的子空间：$0+0=0$，$17\cdot0=0$，封闭性显然成立 |

**R^3 的全部子空间**

| 子空间 | 直观维数 | 说明 |
|---|---|---|
| 整个 $\mathbb{R}^{3}$ | 3 | 空间本身 |
| 过原点的平面 | 2 | 平面必须过原点 |
| 过原点的直线 | 1 | 直线必须过原点 |
| $\{0\}$，即 $\begin{bmatrix}0\\0\\0\end{bmatrix}$ 单独成集 | 0 | 最小子空间 |

> **易错点**：$\mathbb{R}^{2}$ 中过原点的直线虽然「看起来像一维」，但它**不等于 $\mathbb{R}^{1}$**——线中的向量仍有两个分量，而 $\mathbb{R}^{1}$ 的向量只有一个分量。子空间是「集合层面」的概念，不能凭画得像就划等号。

---

## 5. 从矩阵产生子空间：列空间

### 5.1 定义

把矩阵 $A$ 的各列看作向量，它们**全部线性组合**的集合，称为 $A$ 的**列空间（column space）**，记作 $C(A)$。

**为什么线性组合的全体一定是子空间**：两个线性组合相加仍是线性组合；一个线性组合乘上标量仍是线性组合；乘数全取 $0$ 就得到零向量。所以封闭性自动成立。

**算例（黑板）**：取以两列为向量、位于 $\mathbb{R}^{3}$ 中的矩阵

$$A=\begin{bmatrix}1&3\\2&3\\4&1\end{bmatrix},\qquad
c_1=\begin{bmatrix}1\\2\\4\end{bmatrix},\quad c_2=\begin{bmatrix}3\\3\\1\end{bmatrix}$$

$C(A)$ 必须包含：

- $0\cdot c_1+0\cdot c_2=\begin{bmatrix}0\\0\\0\end{bmatrix}$（零向量必须在内）；
- $c_1+c_2=\begin{bmatrix}4\\5\\5\end{bmatrix}$（加法必须在内）；
- $c_1+3c_2$ 等一切 $x_1c_1+x_2c_2$。

**几何图像**：$c_1$ 与 $c_2$ **不在同一条直线上**（不共线），所以它们的组合先填满各自的整条直线，再填满两线之间的所有点——结果是一个**过原点的平面**，即 $\mathbb{R}^{3}$ 的二维子空间。

> **对照**：如果这两个列向量恰好共线，$C(A)$ 就退化成一条**过原点的直线**。列空间到底多大，完全取决于这些列向量彼此「独立」到什么程度。

### 5.2 推广到高维

五个属于 $\mathbb{R}^{10}$ 的向量取全部组合，得到的同样是某个子空间：

- 它不是 $\mathbb{R}^{5}$——这些向量本身有 $10$ 个分量，空间在 $\mathbb{R}^{10}$ 里；
- 若这五个向量彼此独立，得到的是一个过原点的 $5$ 维「平坦」子空间；
- 若它们之间有冗余（比如全在同一条直线上），维度就相应降低。

**统一观念**：手上只有有限几个向量时，不要满足于「几个向量」，而要取它们**全部线性组合**，从而得到一个**空间**。

---

## 6. 下一讲的引子

本讲正式引出的矩阵子空间只有**列空间**。下一讲会把 $Ax=b$ 翻译成子空间语言：

$$Ax=b \text{ 有解}\quad\Longleftrightarrow\quad b \in C(A)$$

即右端项 $b$ 必须落在 $A$ 的列空间里。另一些矩阵子空间（例如**零空间 nullspace**）留待后续讲次展开。

---

## 总结

1. **消除行交换的障碍后才有一切**：一般消元公式是 $PA=LU$；$A=LU$ 只是 $P=I$ 的特例。
2. **置换矩阵三件事**：$P$ 是重排行的单位矩阵；$n\times n$ 共 $n!$ 个；$P^{-1}=P^{T}$（因为 $P^{T}P=I$）。
3. **转置规则与求逆同构**：$(AB)^{T}=B^{T}A^{T}$，$(A^{T})^{T}=A$；对称矩阵即 $S^{T}=S$，识别方式是沿主对角线镜像相同。
4. **$R^{T}R$ 永远对称**：由 $(R^{T}R)^{T}=R^{T}R$ 两行证完；黑板数字算例给出 $R^{T}R=\begin{bmatrix}10&11&7\\11&13&11\\7&11&17\end{bmatrix}$。
5. **向量空间与子空间的判据**：对加法与数乘封闭（等价于对任意线性组合封闭）；由此立刻推出**必须含零向量**——这是最快的排除法。第一象限对加法封闭但对数乘不封闭，故不是空间。
6. **子空间清单**：$\mathbb{R}^{2}$ 只有三类——整个平面、过原点的直线、$\{0\}$；$\mathbb{R}^{3}$ 四类——整个空间、过原点的平面、过原点的直线、$\{0\}$。
7. **列空间的构造方式**：取矩阵的列，再取全部线性组合；黑板两列不共线，因此张成 $\mathbb{R}^{3}$ 中一个过原点的平面；若两列共线则退化为直线。

---

## 参考资料

- [MIT OCW 18.06 Linear Algebra 课程主页（视频、讲义、习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/)
- [本讲页面：Lecture 5: Transposes, permutations, spaces R^n（视频、讲义与习题）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-5-transposes-permutations-spaces-r-n/)
- [下一讲页面：Lecture 6: Column space and nullspace（列空间与零空间）](https://ocw.mit.edu/courses/18-06-linear-algebra-spring-2010/resources/lecture-6-column-space-and-nullspace/)
- [MIT 18.06 官方课程网站（web.mit.edu/18.06）：按主题索引的视频与讲义](https://web.mit.edu/18.06/www/)
- [MIT 18.06 课程要点摘要（mitmath/1806 仓库，含 $PA=LU$ 与子空间要点）](https://github.com/mitmath/1806)
- [Gilbert Strang《Introduction to Linear Algebra》教材官网（目录与配套材料）](https://math.mit.edu/~gs/linearalgebra/)
- [MIT 18.06 线性代数中文笔记（ApacheCN 翻译整理，供中文术语对照）](https://github.com/apachecn/mit-18.06-linalg-notes)

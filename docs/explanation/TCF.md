# TCF Explanation

| 标题翻译 |                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 作者     | Pengcheng Shi; Shaocheng Yan; Yilin Xiao; Xinyi Liu; Yongjun Zhang; Jiayuan Li                                                                                                                                                                                                                                                                                                                                                                 |
| 出版年份 | 12/2024                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| 期刊     | IEEE Robotics and Automation Letters                                                                                                                                                                                                                                                                                                                                                                                                           |
| 期刊等级 | <span style="color: rgb(111, 69, 135);"><span style="background-color: rgb(232, 222, 238);">ㅤㅤ ㅤㅤIF 5.3 ㅤㅤ ㅤㅤ</span></span> ㅤㅤ ㅤㅤ<span style="color: rgb(204, 31, 0);"><span style="background-color: rgb(255, 226, 221);"> ㅤㅤ ㅤㅤSCI Q1 ㅤㅤ ㅤㅤ</span></span> ㅤㅤ ㅤㅤ<span style="color: rgb(111, 69, 135);"><span style="background-color: rgb(232, 222, 238);"> ㅤㅤ ㅤㅤSCI升级版 计算机科学2区 ㅤㅤ ㅤㅤ</span></span> |

***

![TCF algo 1](./assets/TCF_algo_1.png)

第一个RANSAC很简单， model就是满足式(3)的点集，式(3)就是在约束尺度变化量，找一对锚correspondence，然后比较各自与其他点的长度，足够小的就放进consensus set，注意这个式子是通过三角形边长不等式推导得到的。作者并没有说$\tau$是怎么算出来的，尽管他说可以automatically adjust，需要看代码。

![TCF algo 2 page 1](./assets/TCF_algo_2_p1.png)

![TCF algo 2 page 2](./assets/TCF_algo_2_p2.png)

第二个RANSAC的model是在第一consensus set的基础上找满足式(6)-式(8)的点集，首先找俩对锚correspondence，通过角度约束得到第二个consensus set，需要注意的在最坏情况下（$p_m$和$q_m$重合，就是第二张图所示），a最多可以为$2\tau$,这是因为此时$p_m$刚好在$p_iq_i$的中垂线上，c也同理地小于等于$2\tau$；而$\beta$的分子为$\tau$是因为这是对于单个correspondence而言的最大角度偏差。

第三个RANSAC就是Kabsch，没啥好讲的，需要注意的是，他可以直接用在algorithm 3的初始化（下图式1）。

![TCF algo 3](./assets/TCF_algo_3.png)

step1：
通过Kabsch得到初始变换，然后计算得到各correspondence的residual，选个最大的当作初始scale，后文会提到这个叫Cauchy分布的尺度参数。intialize里的w\_1是1，这代表是L2-norm，即所有点一视同仁，无权重。

step3&4：
通过IRLS得到带权重的最小二乘解，具体原理看：<https://sepwww.stanford.edu/data/media/public/docs/sep115/jun1/paper_html/node2.html>，根据得到的解算出估计的最小residual，丢进Cauchy分布里面，如果residual在$3\gamma$（这个“3”可用于估计器的一致性判断，这里应该是个经验值）范围里面（$\gamma$正是Cauchy分布的尺度参数，见下图），那就认为是内点。

![Cauchy Loss](./assets/cauchy_loss.jpg)

step5&6：
选这个内点residual，通过Cauchy Loss更新权重，注意big residual mostly outlier, small residual mostly inliner，这句话你可以看Cauchy Loss的表达式自己悟出来。注意Cauchy Loss会比Huber Loss对big outlier更不敏感，从而更好剔除，具体可以看看step3&4提到的网址。

step7：
scale通过除以衰减因子$\mu$一步一步缩小，直到找到最佳变换，有点像吃鸡里面的“缩毒圈”:>

step8：
进行收敛判断：两个条件满足一个就结束收敛：（1）加权residual平方差小于$\mathbf{e}_{\text{min}}$，其中e\_j就是上一次的估计residual。这个$\mathbf{e}_{\text{min}}$按AI的说法是最小能量阈值，我的理解就是residual的变化量小到一定程度以至于梯度最小，那么就是收敛了。（2）scale比1小。从上图来看，scale到1的内点范围已经很小了，再小就退化成冲激函数（数学里叫Kronecker函数）了！

***

个人观点：感觉可以优化的点在于Kabsch、衰减因子$\mu$和最小能量阈值$\mathbf{e}_{\text{min}}$

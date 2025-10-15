# RGS2

We present new *randomized Garm-Schmidt (RGS)* algorithms inspired by random sketching and reorthogonalizing. While the existing RGS algorithm developed by [Oleg and Grigori, 2022] provides a set of sketched orthonormal basis vectors, our proposed methods can generate $l_2$ orthonormal vectors. Compared to other *classical/modified Gram-Schmidt (CGS/MGS)* algorithms, we can obtain benefits in reduced computational complexity as well as numerical stability. Round-off error analyses are proved in [Randomized Orthogonalization Process With Reorthogonalization](https://onlinelibrary.wiley.com/doi/10.1002/nla.70029).

We give two types of RGS2 algorithms, called **RGS2C** and **RGS2M**, according to reorthogonal process, respectively. To validate the advantage of methods, we carry out some numerical experiments. First, we compare the numerical performance with various GS algorithms. Furthermore, we apply our algorithms to GMRES and GMRES-DR for solving lienar systems. Consequently, the proposed methods exhibit improved numerical performance. 

If you have any enquiry, please contact **Yongseok Jang (yongseok.jang@chonnam.ac.kr)**.

# References

Main references used in ESTstack. See per-algorithm explanation pages for detailed equation-level citations.

---

## Kalman Filter

```bibtex
@misc{wang2023bayesian,
  title        = {贝叶斯滤波重制版},
  author       = {忠厚老实的老王},
  howpublished = {Bilibili video series},
  year         = {2023},
  url          = {https://www.bilibili.com/video/BV1Fj41157Ey/},
  note         = {Bayesian filtering, Kalman filter, Markov processes,
                  nonlinear KF (EKF/UKF), particle filter, and resampling}
}

@misc{zhang2024sdm366,
  title        = {{SDM366} -- Optimal Control and Estimation},
  author       = {Zhang, Wei},
  institution  = {Southern University of Science and Technology (SUSTech)},
  howpublished = {Course, Spring 2024},
  year         = {2024},
  url          = {https://github.com/clearlab-sustech/SDM366-Sp24},
  note         = {State-space modeling, least-square estimation, system identification,
                  LQR, observer design, and Kalman filter}
}

% Graduate-level reference; requires VERY SOLID linear algebra, probability (including stochastic
% process) and Lie group background. Highly recommended once fundamentals are REALLY in place.
@book{barfoot2024state,
  title     = {State Estimation for Robotics},
  author    = {Barfoot, Timothy D.},
  edition   = {2},
  publisher = {Cambridge University Press},
  year      = {2024},
  isbn      = {978-7-5693-0624-8},
  url       = {http://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser24.pdf},
  note      = {Bayesian estimation, Kalman filter, EKF, UKF, batch optimization,
               Lie groups for 3D estimation, Gaussian processes on trajectories}
}
```

## ESKFOM

```bibtex
@article{sola2017quaternion,
  title   = {Quaternion kinematics for the error-state {K}alman filter},
  author  = {Sol{\`a}, Joan},
  journal = {arXiv preprint arXiv:1711.02508},
  year    = {2017},
  url     = {https://arxiv.org/abs/1711.02508},
  note    = {Local/global perturbation ESKF with Hamilton quaternions;
             IMU error equations and low-cost IMU noise modeling}
}

@article{sola2018micro,
  title   = {A micro {L}ie theory for state estimation in robotics},
  author  = {Sol{\`a}, Joan and Deray, Jeremie and Atchuthan, Dinesh},
  journal = {arXiv preprint arXiv:1812.01537},
  year    = {2018},
  url     = {https://arxiv.org/abs/1812.01537},
  note    = {Lie group operations (exp, log, Adjoint, Jr/Jl) on
             SO(3), SE(3), SE_2(3) and Bundle groups}
}

@techreport{blanco2010tutorial,
  title       = {A tutorial on {SE(3)} transformation parameterizations
                 and on-manifold optimization},
  author      = {Blanco-Claraco, Jos{\'e} Luis},
  institution = {University of Malaga},
  year        = {2010},
  url         = {https://ingmec.ual.es/~jlblanco/papers/jlblanco2010geometry3D_techrep.pdf},
  note        = {SE(3) representations, Jacobians, and manifold-based optimization}
}

@article{qiu2019monocular,
  title   = {Monocular Visual-Inertial Odometry with an Unbiased Linear
             System Model and Robust Feature Tracking Front-End},
  author  = {Qiu, Xiaochen and Zhang, Hai and Fu, Wenxing and Zhao, Chenxu and Jin, Yanqiong},
  journal = {Sensors},
  volume  = {19},
  number  = {8},
  pages   = {1941},
  year    = {2019},
  note    = {Unbiased linear system model for VIO}
}

@misc{sola2021lie_video,
  title        = {Lie theory for the roboticist},
  author       = {Sol{\`a}, Joan},
  howpublished = {YouTube, Institut de Rob\`otica i Inform\`atica Industrial, CSIC-UPC},
  year         = {2021},
  url          = {https://www.youtube.com/watch?v=nHOcoIyJj2o},
  note         = {Visual introduction to Lie groups, Lie algebra, exp/log maps,
                  Adjoint, and left/right Jacobians for robotics ESKF}
}
```

---

## CT (Coordinated Turn) Exact Modeling

```bibtex
@misc{hammarstrand_ssy345,
  title        = {{SSY345} -- Sensor Fusion and Nonlinear Filtering},
  author       = {Hammarstrand, Lars},
  institution  = {Chalmers University of Technology},
  howpublished = {edX online courses},
  note         = {Nonlinear motion models (CV, CA, CT, CTRV),
                  EKF, UKF, particle filter, and IMM}
}
```

---

## BTC Descriptor 

```bibtex
@article{yuan2024btc,
  title   = {{BTC}: A Binary and Triangle Combined Descriptor
             for 3-{D} Place Recognition},
  author  = {Yuan, Chongjian and Lin, Jiarong and Liu, Zheng
             and Wei, Hairuo and Hong, Xiaoping and Zhang, Fu},
  journal = {IEEE Transactions on Robotics},
  volume  = {40},
  pages   = {1580--1599},
  year    = {2024},
  doi     = {10.1109/TRO.2024.3353076},
  note    = {Binary triangle descriptor for LiDAR place recognition;
             geometric consistency voting for robust matching}
}
```

---


### Two-Stage Consensus Filtering

```bibtex
@article{shi2024ransac,
  title   = {{RANSAC} Back to {SOTA}: A Two-Stage Consensus Filtering
             for Real-Time {3D} Registration},
  author  = {Shi, Pengcheng and Yan, Shaocheng and Xiao, Yilin
             and Liu, Xinyi and Zhang, Yongjun and Li, Jiayuan},
  journal = {IEEE Robotics and Automation Letters},
  volume  = {9},
  number  = {12},
  pages   = {11881--11888},
  year    = {2024},
  doi     = {10.1109/LRA.2024.3502056},
  note    = {Coarse compatibility graph voting + fine RANSAC verification
             for real-time 3D registration}
}
```

---

## Code Refinement

```bibtex
@article{deray2020manif,
  title   = {manif: A small header-only library for {L}ie theory},
  author  = {Deray, Jeremie and Sol{\`a}, Joan},
  journal = {Journal of Open Source Software},
  volume  = {5},
  number  = {46},
  pages   = {1371},
  year    = {2020},
  url     = {https://github.com/artivis/manif},
  note    = {C++ Lie groups: SO2, SO3, SE2, SE3, SE_2_3, Bundle, Rn;
             automatic Jacobians via plus/minus}
}

@misc{deray2021kalmanif,
  title   = {kalmanif: A small collection of {K}alman Filters on {L}ie groups},
  author  = {Deray, Jeremie},
  year    = {2021},
  url     = {https://github.com/artivis/kalmanif},
  note    = {EKF, IEKF, UKFM, RTS smoother on manifolds;
             SE_2(3) IMU strap-down ESKF reference (demo\_se\_2\_3.cpp)}
}

@book{gao2024slam,
  title     = {自动驾驶与机器人中的{SLAM}技术},
  author    = {高翔},
  publisher = {电子工业出版社},
  year      = {2024},
  isbn      = {978-7-121-45878-1},
  note      = {Covers LiDAR/visual SLAM, Kalman filtering, graph optimization,
               and point cloud registration with modern C++ implementation}
}
```

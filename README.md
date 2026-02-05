This repository provides the complete numerical implementation supporting the results and figures reported in the accompanying manuscript. Its primary purpose is to ensure transparent, verifiable, and long-term reproducibility of the published analysis.

The code implements a physics-based model of bow mechanics, capturing the coupled effects of limb geometry, rotational stiffness, string kinematics, and energy transfer for three distinct bow architectures: a self-bow, a composite recurve bow, and an asymmetric Japanese Yumi. All force–draw relations, energy calculations, efficiency measures, and sensitivity analyses reported in the manuscript are generated directly from this implementation.

The repository is intentionally structured so that all figures can be reproduced with a single command, without the need for external data files or manual parameter tuning. Numerical solvers, integration schemes, and random sampling procedures are fully specified and deterministic, enabling exact replication of results across systems.

This code is provided primarily for reviewers and researchers interested in verifying the analysis, reproducing the figures, or extending the model to related historical or biomechanical bow configurations. Familiarity with Python is sufficient to run the simulations; detailed knowledge of the internal implementation is not required for reproduction.

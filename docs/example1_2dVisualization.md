# Example 1: 2D visualization

- spatial(2D) & UMAP of one sample

  In this example, `PANEL 1` and `PANEL 2` are synchronized by `feature` and `sample`, which means their sample and feature to plot will be kept the same. Then you can check different genes without adjust them separately.

  <img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20250103162409550.png" alt="image-20250103162409550" style="zoom:50%;" />

  

- spatial(2D) of different samples

  In this example, `PANEL 1` and `PANEL 2` are synchronized by `feature`, `highlighting` and `embedding`, but no `sample`. Therefore, we could check, in multiple samples simultaneously, the items in metadata or specific genes on same embedding. These samples could be from replicates, different stage, or different datasets.

  <img src="https://pic-md-1259550128.cos.ap-nanjing.myqcloud.com/image-20250110182304449.png" alt="image-20250110182304449" style="zoom:50%;" />


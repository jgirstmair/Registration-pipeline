This project tries to perform batch registration of embryos.

<img src="https://user-images.githubusercontent.com/70565347/126686183-c9ce7a43-bb71-4208-a939-2758eccc998e.png" width=400>


### TODOs

- [x] Perform isotropic resolution of sample
- [x] Read Pixel Size form Image in Python
- [x] Fix the 16 bit auto conversion in Isotropic Sampling in Python (Keep 8 bit)
- [x] Perform 3D pre-rotation in Fiji
- [x] Perform 3D pre-rotation in Python
- [x] Do not save isotropic as RGB
- [x] Apply modification to bring to canonical orientation
- [x] Set the final size such that the warped image is centred on the specimen
- [x] Turn head to a consistent side (done with the intensity of *phalloidin* detections)
- [x] For pre-registration, write a batch script to process all files
- [ ] Use detections only once (currently, for the head consistent alignment, we detect them again)
- [ ] Extend current python script to batch-process rigid and FFD registration
- [ ] Establish Averaging Pipeline
- [ ] Not all `Macrostomum` are perfectly aligned (~ 5/15: were aligned vertically with random head position)
- [ ] Investigate some outlier filtering for detections before doing PCA
- [ ] Produce results on larva (more spherical with lobes)
- [ ] Perform alignment using keypoint detection (pose estimation) on larvae (for example, see [this](https://www.nature.com/articles/s41593-018-0209-y))
- [ ] Perform Noise2Void on all noisy images
- [x] Grow Images by 1.3 (currently set to 1.2) 
- [ ] Rotate randomly and check if you can recover the correct alignment. First test failed.
- [ ] PCA rotation matrix and the pre-registration should be applied on the original (not down-sampled) images.
- [ ] Work on Bruno's data. 
- [ ] Register all `phalloidin` to some `Phalloidin` source image for `Echinoplana` 
- [ ] Investigate if the process of rotation generates empty background which confuses PCA (for skewed embryo)
- [ ] Compare denoised images from cluster with the ones generated before ( Send images to Johannes )
- [ ] Create first test of registering allneuropeptides againsta sensible source image.
- [ ] See if normalization prior to PCA makes it more robust to outliers.
### Learnings

- Preregistration works betteron reference channel 1 (`phalloidin`). With `DAPI`, we have more variation and flipped left-right heads are more often wrong. 

- First vector of PCA does correspond to longest axis or axis of largest variance, but it is senstive to outliers. We noticed that with `Macrostomum`, the `DAPI` channel didn't perfrom so well because of the outliers - hence, we switched to `Phalloidin` for detections and PCA.

- To get a consistent detection ogf head, we tried some algorithms - but what worked the best finally was using the intensity at detections. We use the center plane at the specimen center and compare the intensities on the left and right (this shoudl correspond to the AP-axis) and should gve a bias along all bilaterians. Note that for `echinoplana`, the pharynx lies closer to the trunk while for the `macrostomum`, it is closer to the head, hence the algorithm is flipped but consistent. Also, it seems to work with DAPI in `echnioplana`.

- To apply the transformation matrix from PCA, we ralised that the rotation should occur around the embryo centre. [This](https://simpleitk.readthedocs.io/en/v1.2.4/Documentation/docs/source/fundamentalConcepts.html) page explains how the transformation is done 

- We expect the rigid registration to work perfect. 

- We also noticed that using the nervous system channel performs 100 % correctly for `echninoplana` and performs at 12/16 on `macrostomum`. The fact that the nervous system shows higher intensity in the brain can be used to align the head consistently. 

#### Important files

- python/Registration_Pipeline_macrostomum.ipynb
- python/pythonScripts/preregisterAnimals.py


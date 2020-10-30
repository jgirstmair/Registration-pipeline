This project tries to perform batch registration of embryos.


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
- [ ] 

### Learnings

- First vector of PCA does correpsond to longest axis or axis of largest variance, but it is senstive to outliers. We noticed that with `Macrostomum`, the `DAPI` channel didn't perfrom so well because of the outliers - hence, we switched to `Phalloidin` for detections and PCA.

- To get a consistent detection ogf head, we tried some algorithms - but what worked the best finally was using the intensity at detections. We use the center plane at the specimen center and compare the intensities on the left and right (this shoudl correspond to the AP-axis) and should gve a bias along all bilaterians. Note that for `echinoplana`, the pharynx lies closer to the trunk while for the `macrostomum`, it is closer to the head, hence the algorithm is flipped but consistent. Also, it seems to work with DAPI in `echnioplana`.

- 





# Gaussian Mixture Midway-Merge for Object SLAM with Pose Ambiguity

**Authors : [Jae Hyung Jung](https://sites.google.com/view/lastflowers), and Chan Gook Park**

### 1. Overview

This is a MATLAB script that implements the Gaussian-mixture invariant EKF with midway-merge in the Chapter 6 of the paper.  This demonstration reproduces results of the 0022 sequence in the YCB-Video dataset. Pose detection estimates from [CosyPose](https://github.com/ylabbe/cosypose) with a single view are prepared in the `cosypose/0022/` folder with the name format as `obj_id_imageindex.npy`. Note that we use pretrained weights from the authors to estimate camera-object relative poses.

[![Video Label](http://img.youtube.com/vi/EvtW-mb8YK8/0.jpg)](https://youtu.be/EvtW-mb8YK8)



### 2. Run (YCB-Video example)

* Download the YCB-Video toolbox: https://github.com/yuxng/YCB_Video_toolbox
* Download the YCB-Video dataset 0022 sequence for a quick demonstration.
* Install npy-matlab to read npy files from the pose detector: https://github.com/kwikteam/npy-matlab
* Set your path in the `main.m` and run.

```matlab
path_toolbox = 'path_to_ycb_video_toolbox';
path_dataset = 'path_to_dataset';
```

* To visualize the estimated masks with the current image, please set `isPlot = 1;`. But, note that this will slow the process due to the drawing.



### 3. Citation

If you feel this work helpful to your academic research, we kindly ask you to cite our paper :

```
@article{Midway_RAL,
  title={Gaussian Mixture Midway-Merge for Object SLAM with Pose Ambiguity},
  author={Jung, Jae Hyung and Park, Chan Gook},
  journal={IEEE Robotics and Automation Letters},
  year={2022},
  publisher={IEEE}
}
```



### 4. Acknowledgements

This research was supported by the National Research Foundation of Korea funded by the Ministry of Science and ICT, the Republic of Korea. (NRF-2022R1A2C2012166)



### 5. License

Our source code is released under MIT license. If there are any issues in our source code please contact to the author (lastflowers@snu.ac.kr).




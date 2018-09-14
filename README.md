# GMSimViz: Automated 3D Visualization of Ground Motion Simulation with Generic Mapping Tools (GMT)

## Overview
GMSimViz is an automation tool that produces an animated 3D visualization of geological faults, ground motion and other earthquake related data. It uses the Generic Mapping Tools(GMT) to create individual frames, and join them by FFmpeg to create a movie file.

A fully featured result is available at https://youtu.be/qZkOTI4x_cc

<table style="width:%">
  <caption>Screenshots of a GMSimViz-generated animation</caption>
  <tr>
    <th><img src="figures/kaikoura_fault_slip_distribution.jpg" width="400"></th>
    <th><img src="figures/kaikoura_ground_motion_hitting_wellington.jpg" width="400"></th>
  </tr>
  <tr>
    <td><img src="figures/kaikoura_peak_ground_velocity.jpg" width="400"></td>
    <td><img src="figures/kaikoura_landslide_susceptibility.jpg" width="400"></td>
  </tr>
  <tr>
    <td><img src="figures/kaikoura_liquefaction_susceptibility.jpg" width="400"></td>
    <td><img src="figures/kaikoura_transport_recovery.jpg" width="400"></td>
  </tr>
</table>



## Installation

### Dependencies
The GMSimViz depends on the following software packages.

* Python (>=2.7) tested with 2.7 and 3.6. Which version do I have? Run `/usr/bin/env python --version`.
* Numpy (https://www.numpy.org https://pypi.org/project/numpy/) for Python version (>= 1.10) tested with 1.10, 1.13 and 1.14
* MPI4Py (https://bitbucket.org/mpi4py/ https://pypi.org/project/mpi4py/) for Python version and mpi backend (https://www.open-mpi.org dependency of mpi4py)
* H5Py (http://www.h5py.org/ https://pypi.org/project/h5py/) for Python version and hdf5 backend (https://support.hdfgroup.org/HDF5/ dependency of h5py)
* GMT (>=r19922) requires release after 5.4.3 (currently unavailable) versions prior to r19922 will have bugs but GMSimViz is designed to work with GMT (>=5.2) tested with Ghostscript 9.18, 9.21 was found to produce glitches
* FFMpeg built with image2/png, h.264 encoder support (standard installation) tested with version 3.3
* gawk

### Linux
Install the packages for the above mentioned dependencies. Names vary between distribution repositories.  On Ubuntu or other Debian-based Linux distributions, most of the packages can be installed by the following command. If your system default Python is version 3, install the `python3-*` packages instead of `python-*`.
```shell
sudo apt install python2.7 python-numpy python-mpi4py python-h5py ffmpeg gawk
```
If you are using a Debian based distribution, you will also need to have installed the separated -dev packages in order to compile GMT.
```shell
sudo apt install libpng-dev
```
#### GMT

GMT version earlier than 5.2 has some outstanding issues to work reliably with GMSimViz. A version later than 5.4.3 is recommended, and at the time of writing, the latest development source obtained directly from the  repository (A release r19922 fully verified).  Instructions for GMT installation is available at http://gmt.soest.hawaii.edu/projects/gmt/wiki/BuildingGMT . You should compile with PNG support enabled (check configure output, make sure headers are available).

If you want to use the latest stable release at time of writing, version 5.4.3 will draw outlines of the fault plane and hypocentre location on the surface rather than at depth (bug).

Make sure the `gmt` binary is in the PATH. If `gmt` is installed in `/usr/local/bin` and `gmt` is not already available in the PATH, add it to the PATH:
```shell
export PATH=$PATH:/usr/local/bin
```
You can add the line above to `~/.bashrc` to make it persistent.


#### GMSimViz

Download the plotting resource GMSimViz_resources.zip file from https://goo.gl/mYFCQn
Extract this file to where GMSimViz is located so that GMSimViz/resources exists.

Install with the following command.
```shell
sudo python setup.py install
```
In case of no root priviledge, `--user` option can be used.

```shell
python setup.py install --user
```
You may use `python2`, `python3`, `python2.7` etc. instead of `python` where you want to use a specific verision.

### Other systems
Currently unsupported (we cannot give instructions but if the dependencies are met, GMSimViz should work).


## Sample Data
Sample data is provided, it contains smaller data sets from multiple earthquakes to make a feature-full demo. 

The most simplest output is a single frame (image) facing the fault plane. This case only requires the SRF file as the parameter:
<img src="figures/fault_perspective.png" width="400">
```shell
gmsimviz sample_data/fault.srf
```

An animation can be created from an SRF that shows the rupture propagating. It can be expected to take around 12 minutes to complete on a decent personal machine. Here is an example using 7 slave processes for a machine with 8 cores (including hyperthreading). An example of what is expected is available at https://youtu.be/nTtq_DxGVUM
```shell
gmsimviz sample_data/fault.srf -a --crude -n7 --title "Fault Animation" --dpi 120 --downscale 1
```

The fully featured sample animation below takes about 45 minutes with -n7 on a v4 E5-2620 (2.1 GHz). Slave processes are spawned after about 3 minutes and frames will start appearing in the temporary folder (in the same folder as the SRF file, sample_data/_GMT_WD_PERSPECTIVE_<random sequence>). When complete, the video is available in the working directory and the temporary folder is removed by default.
```shell
gmsimviz sample_data/fault.srf -a --crude -n7 --title "Sample Animation" --dpi 120 --downscale 1 -x sample_data/xyts.e3d --liquefaction-s sample_data/liquefaction_s.hdf5 --liquefaction-p sample_data/liquefaction_p.hdf5 --landslide-s sample_data/landslide_s.hdf5 --landslide-p sample_data/landslide_p.hdf5 --paths sample_data/transport
```

Here are details for command options:

- `srf_file` first position, not optional. Path to a fault in standard rupture format containing the plane header (https://scec.usc.edu/scecpedia/Standard_Rupture_Format).
- `-a` create animation instead of static image showing the slip distribution only. Will animate the progression of slip if no other input data are given.
- `--crude` this will produce a quick result. No topography or roads are plotted, low resolution coastlines are used and overlays have lower resolutions. This is required for faults outside of the New Zealand area because the higher resolution topography is only available here. Alternatively, replace the topography files with your own.
- `-n28` number of slave (worker) processes to start. 7 is a good choice for systems with only 8 cores. It is recommended to run on a machine with many cores.
- `-f10` set the framerate of the animation.
- `--title` title on the movie, the default is the basename of the SRF file.
- `--dpi` dpi of output. Frames are 16 inches x 9 inches so 240 will produce 4k output, 120 for FullHD.
- `--downscale` render at higher resolution and then downscale. Prevents jitter in object positioning. 8 is ideal for a smooth result.
- `-x` XYTS file. This is an output created by the EMOD3D software that simulates ground motion. This provides the simulation domain and ground motion data.
- `--liquefaction-s` and `--landslide-s` Liquefaction and/or landslide suceptibility HDF5 file. Must have data under `model` as in given sample data.
- `--liquefaction-p` and `--landslide-p` As above for the probability from given event.
- `--paths` Transport network in GMT path format with optional GMT legend files as given in sample data.
- `--gm-cut` Cut the time series of the xyts file given with `-x` short to this many seconds.
- `-t` Sets the time (s) for major animation transitions such as when originally facing the slip distribution.
- `-m` Sets the minor tranition time (s) for transitioning between different segments of the animation.
- `-p` Sets the time (s) to pause for when showing an image such as fault under ground and liquefaction/landslide images.
- `-d` Sets the delay (s) at the start of the animation before the camera starts moving.
- `-e` Sets end delay (s) at the end of the animation to prevent sudden end after moving the camera.
- `-r` Fix the rotation of the camera to bearing (deg). eg: 0 to face north always.
- `--temp` Continue from the given previous temporary directory, will not re-render frames which are complete.
- `-k` Keeps the temporary directory (containing all frames) once animation is complete.

To list the parameters, run:
```GMSimViz.py -h```

The SRF file format is available here: https://scec.usc.edu/scecpedia/Standard_Rupture_Format. The XYTS file format is described [here](./XYTS.md).

## Thanks

The QuakeCoRE team would like to thank:

* Eric Thompson from USGS for collaborating with us and providing data and software relating liquefaction and landslide analysis.
* The Generic Mapping Tools team, especially Paul Wessel for GMT support and quickly fixing issues we ran into while using GMT.

# Spy(name yet to be selected)

Spy is a 3D to 2D correspondence estimator from Indshine. 
It works together with Indshine's [Image Viewer](http://www.indshine.com)
and [AWS Lambda](https://aws.amazon.com/lambda/)service. It can also be run on a 
local PC.

## Getting Started

These instructions will get you a copy of *Spy* up and running on 
your local or cloud machine for development and testing purposes.

### Prerequisites

Spy is alomost a standalone software and only depends on few basic libraries
1. [numpy](http://www.numpy.org/) is used to convert complex calculations into
simple matrix multiplication task.
2. [utm](https://pypi.org/project/utm/) is to easily transform geographic 
coordinate system into projected coordinate system and vice versa.
3. [math](https://docs.python.org/2/library/math.html) is used to do simple 
mathematical operations.

Spy also need a image location file containing information like Lat, Long, 
Elevation, Omega, Phi, Kappa or Rotational Matrix parameters.
Output can be obtained from [openMVG](https://github.com/openMVG/openMVG), 
[openSFM](https://github.com/mapillary/OpenSfM) or 
[Agisoft](http://www.agisoft.com/). 
* Currently openMVG and oepnSFM format are not supported.


### Installing

A step by step series of examples that tell you how to get a development 
env running on your local PC

1. Simply clone this repository. To get the latest version of Spy, work with 
**develop** Branch while **master** branch gives you earlier but stable version.

2. There are three function that come along with Spy. These are
**Matrix_Projection** , **BackProject** , **Lambda_function**

3. **Matrix_Projection** is used to project 100s or 1000s of 3D points to a single 
image at once. These 3D points must be in 2D format. It is optimised to do 
vectorised operations of large number of points. Output will be pixel X and 
pixel Y location of 3D points.  

4. **BackProject** takes single point as input and project that point to all the 
available images. It gives pixel X and pixel Y location in all the visible 
images. If the X and Y location is greater than column and row of image then it
is automatically discarded and marked as not visible in that image.

To run Matrix_Projection, simply 
1.  create 2D matrix containing information of lat , long , elevation of points 
to be projected 

2.  String varible of location of camera file containing latitude, longitude , 
elevation , and rotational matrix information. If omega, Phi, Kappa is 
available then use eulerAnglesToRotationMatrix_wiki function to 
create a rotational matrix.

3. String variable containing name of image to which all 3D points are to 
be projected

then use import Matrix_Projection funtion to import in environment 

```
Example -  import Matrix_Projection as mp
           Pixels = mp.BackProject(lat,long,elevation,'camera.txt','DSC_0323.JPG')
```

Output will contain 2 dimension matrix, first element will give x coordinates 
and second contains y coordinates. Each x and y coordinates is itself a 2D matrix
of same size as of Lat, Long provide. x[1] will coorespondes to x coordinate with 
respect to lat[1], long[1] and elevation[1].

## Running the tests

To run test data use **camera.txt** file given in the dataset. Run in python 2.7 
as follows:

```
long =  [[77.1765467731,77.2065467731],[77.1265467731,77.2465467731]];
lat = [[31.12878750,31.438878750],[31.298878750,31.198878750]];
elevation = [[2107.37,2067.37],[2137.37,2167.37]];
camera_file = 'camera.txt';;
camera_name = 'DSC04045_geotag';

import Matrix_BackProject as mp
coordinates = mp.BackProject(lat,long,elevation,camera_file,camera_name);
```
One can also check [sampe](https://docs.google.com/spreadsheets/d/12dXUMDaMuzXPH__3mwHE2zCIMFH0Rdx3lNIwuJUOeww/edit?ouid=111414247711934633785&usp=sheets_home&ths=true)
document which contains this calculations and verify its results.

### Break down into end to end tests

This test data has project files Shimla region. Total 49 photos were captured 
in mountaneous region. Agisoft is used for processing of these images and 
camera.txt is exported from it. 
Camera.txt contains lat,long,ele,omega,phi,kappa,RotationalMatrix(r11,r12,...)

## Deployment

Clone this repository and follow running test example

## Built With

* [Python 2.7](https://www.python.org/download/releases/2.7/) - The language used

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Manish Sahu** - Product Development

## License

This project is not licensed yet




Dependencies 
1. Numpy 
2. UTM
3. Math

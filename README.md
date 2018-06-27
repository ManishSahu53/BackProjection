This python 2.7 script will convert 3D point to 2D image plane.

Project consist of Lambda function which is executed with aws Lambda call.

Inputs will be
1. X,Y - coordinate of 3D point in world coordinate system. This might be in 
projected or in Geographic coordiante sytem.

2. Camera file - This file contains information about camera location, 
    orientation and rotational matrix. This file is obtained by exporting
    cameras in Agisfot photoscan
    Format of the file is as follows 
PhotoID, X, Y, Z, Omega, Phi, Kappa, r11, r12, r13, r21, r22, r23, r31, r32, r33

Where photo Id will contain name of camera/photo, 

X,Y,Z are coordinate of camera centre. 

Omega, Phi, Kappa is orientation of that camera


r11,r12.. r33 are elements of 3x3 rotational matrix. Here one point is to be 
noted that Rotational Matrix will be 


[ r11 ,  r12 ,  r13]


[-r21 , -r22 , -r23]


[-r31 , -r32 , -r33]



Dependencies 
1. Numpy 
2. UTM
3. Math

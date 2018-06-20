
def lambda_handler(event,context):
    x = event.get('x')
    y = event.get('y')
    z = event.get('z')
    camera_file_location = event.get('camera_file_location')

    def power(x,p):
        return math.pow(x, p)
    
    def get_txt_file(filename):
       f = open(filename).read()
       matrix = [item.split() for item in f.split('\n')[:-1]];
       return matrix
    
    def create_rot_transl_matrix(matrix):
        rot = np.array(matrix[8:17]); # For agisoft 8:17
        rot_mat = np.resize(rot,(3,3));
        t = np.array(matrix[2:5]);
        t_mat = np.resize(t,(3,1));
        str_mat=np.hstack((rot_mat,t_mat));
        mat = str_mat.astype(np.float)
        return mat   
        
    def create_calibration_matrix(f,cx,cy):
        calibration_matrix = np.zeros((3,3));
        calibration_matrix[0,0] = f;
        calibration_matrix[1,1] = f;
        calibration_matrix[0,2] = cx;
        calibration_matrix[1,2] = cy;
        calibration_matrix[2,2] = 1;
        return calibration_matrix
    
    def eulerAnglesToRotationMatrix_wiki(theta):
        R_x = np.array([[1, 0                 ,                   0], 
                        [0, math.cos(theta[0]), math.sin(theta[0])],
                        [0, -math.sin(theta[0]), math.cos(theta[0])]
                        ])
                      
        R_y = np.array([[math.cos(theta[1]),    0,      -math.sin(theta[1])  ],
                        [0,                     1,      0                   ],
                        [math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                        ])
                     
        R_z = np.array([[math.cos(theta[2]),    math.sin(theta[2]), 0],
                        [-math.sin(theta[2]),    math.cos(theta[2]),  0],
                         [0,                                    0,   1]
                        ])
                                         
        R = np.dot(R_z, np.dot( R_y, R_x ));
        return R  
        
    def create_undistort_matrix(x,y,cx,cy,focal,k1,k2,k3,p1,p2):
        x=(x-cx)/focal;
        y=(y-cy)/focal;
        r = math.sqrt(x*x+y*y);
        x_ud = x*(1+ k1*power(r,2) + k2*power(r,4) + k3*power(r,6)) + 2*p1*x*y + p2*(power(r,2) + 2*power(x,2));
        y_ud = y*(1+ k1*power(r,2) + k2*power(r,4) + k3*power(r,6)) + 2*p2*x*y + p1*(power(r,2) + 2*power(y,2));
        return (x_ud*focal)+cx,(y_ud*focal) +cy   
      
    import numpy as np
    import utm
    import math 
    
    camera_type = "sonya6000"
    #data_loc = "/home/indshine-2/Downloads/Link to SfM/BackTracing/"; 
    camera_file = camera_file_location; #camera Location File
    no_of_markers = 1; # One point at a time;
    camera_data = get_txt_file(camera_file);
    X = np.ones(3);
    C = np.ones(3);
    
    rot_mat = np.zeros((4,4));
    rot_mat[3,3] = 1;
    
    if camera_type == "sonya6000":
        # All units in pixels
        focal = 3994.56722706884;
        k1=-0.0606034595363531;
        k2=0.0472416636682627;
        k3=-0.00970597729150313;
        p1=0.00161210027713506;
        p2=-0.00259257606444914;  
        camera_res = [4000,6000]; #Height and Width
        cx = 15.2451750170639 + camera_res[1]/2;
        cy = -40.4380848998112 + camera_res[0]/2;
        
    elif camera_type =="DJI_P4P":
        focal = 3994.56722706884;
        k1= 0.00298599;
        k2= -0.00769116;
        k3= -0.0079115;
        p1=-0.000129713;
        p2= 0.000221193;  
        camera_res = [3648,5472]; #Height and Width
        cx = 15.0488 + camera_res[1]/2;
        cy = 36.8069 + camera_res[0]/2;
    
    calibration_matrix = create_calibration_matrix(focal,cx,cy);
    no_of_cameras = len(camera_data) -2;
    
    image_name=[];
    x_coord = [];
    y_coord = [];

    calibration_matrix = create_calibration_matrix(focal,cx,cy);
    no_of_cameras = len(camera_data) -2;
    
    for cams in range(no_of_cameras):
        for i in range(1):
           
            rot_ag = np.array(camera_data[cams+2][7:16]); # For agisoft 7:16
            rot_ag = np.resize(rot_ag,(3,3));
            rot_ag = rot_ag.astype(np.double);
            rot_ag = np.transpose(rot_ag); # relative to camera
            

            
            theta =  np.array(camera_data[cams+2][4:7]); 
            theta = theta.astype(np.double);
            theta = 3.14*theta/180; # Degree to Radian
                        
            
            C[0:3] = np.array(camera_data[cams+2][1:4]);
            C = C.astype(np.double);
            C_utm= utm.from_latlon(C[1], C[0]);
            
            C[0] = C_utm[0];
            C[1] = C_utm[1];
            
            X[0] = x;
            X[1] = y;
            X[2] = z;
            
            rel_loc = X-C; # Relative location wrt to camera
            
            projection_matrix_ag = np.dot(rot_ag, rel_loc); # Projection Matrix
            
            pixel_coord_ag = np.dot(calibration_matrix,projection_matrix_ag); # Distorted Pixel Location 
            
            pixel_coord_ag = pixel_coord_ag/pixel_coord_ag[2]; # Back to Homogenous Coordinate system
            
            # Undistorted Pixel Location
            pixel_coord_ag_und = create_undistort_matrix(pixel_coord_ag[0],pixel_coord_ag[1],cx,cy,focal,k1,k2,k3,p1,p2);
            
            if min(pixel_coord_ag_und) > 0 and pixel_coord_ag_und[0] <camera_res[1] and pixel_coord_ag_und[1] <camera_res[0]:
                image_name.append(camera_data[cams+2][0])
                pixel = pixel_coord_ag_und[0:2];
                x_coord.append(pixel[0]);
                y_coord.append(pixel[1]);
                
    data = [{"Image_Name": name, "x_coordinate": x_coordinate, "y_coordinate": y_coordinate} for name, x_coordinate, y_coordinate in  zip(image_name,x_coord,y_coord)] ;
    return data

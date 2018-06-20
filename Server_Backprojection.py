def back_projection(x,y,z,camera_file_location,image_location):

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
    
    def eulerAnglesToRotationMatrix(theta):
        R_x = np.array([[1,         0,                  0                   ],
                        [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                        [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                        ])
                      
        R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                        [0,                     1,      0                   ],
                        [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                        ])
                     
        R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                        [math.sin(theta[2]),    math.cos(theta[2]),     0],
                        [0,                     0,                      1]
                        ])
                                         
        R = np.dot(R_z, np.dot( R_y, R_x ));
        return R   
        
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
        focal = 3894.29210499142;
        cx = -43.1955683312423;
        cy = 18.7344728618194;
    camera_res = [4000,6000]; #Height and Width
    
    calibration_matrix = create_calibration_matrix(focal,cx,cy);
    no_of_cameras = len(camera_data) -2;
    image_name =[]; 
    image_coord=[];
    #for cams in range(no_of_cameras):
    for cams in range(no_of_cameras):
        for i in range(no_of_markers):
            
            #rot = np.array(camera_data[cams+2][7:16]); # For agisoft 7:16
            #rot_mat[0:3,0:3] = np.resize(rot,(3,3));
            #rot_mat = rot_mat.astype(np.double);
            theta =  np.array(camera_data[cams+2][4:7]); 
            theta = theta.astype(np.double);
            theta = 3.14*theta/180; # Degree to Radian
            
            rot = eulerAnglesToRotationMatrix(theta);
            rot = np.transpose(rot); # Camera Rotation to Rotational matrix
            #rot_mat[0:3,0:3] = rot;
            
            C[0:3] = np.array(camera_data[cams+2][1:4]);
            C = C.astype(np.double);
            C_utm= utm.from_latlon(C[1], C[0]);
            
            C[0] = C_utm[0];
            C[1] = C_utm[1];
            X[0] = x;
            X[1] = y;
            X[2] = z;
            
            rel_loc = X-C;
            #rel_loc[3] =1;
            
            projection_matrix = np.dot(rot, rel_loc);
            pixel_coord = np.dot(calibration_matrix,projection_matrix);
            pixel_coord = pixel_coord/pixel_coord[2];
            positive_condition =  pixel_coord>0;
            if positive_condition.all() == True:
                image_name.append(camera_data[cams+2][0])
                image_coord.append(pixel_coord);
                
    data = [{"Image_Name": name, "Corrdinate": coord} for name, coord in zip(image_name, image_coord)] ;
    return data
    
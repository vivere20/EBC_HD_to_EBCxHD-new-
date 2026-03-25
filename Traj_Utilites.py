import numpy as np

def make_envirs(numenv): 
    possible_types = ['Square','Hall','Circle',]
    the_types = possible_types[0:numenv]
    env_types = []
    boundary_points = dict()
    
    area = 10000 #cm^2
    # Square environment: 1 m x 1 m box
    # Center of the environment: approximately 0,0
    # ~1 cm resolution for boundary points
    
    for curenv in the_types:
        print(curenv)
        env_types.append(curenv)
        
        if curenv == 'Square': 
            L = np.sqrt(area)           # side length
            step=1.0
            minc, maxc = -L/2, L/2      # square corners
            coords = np.arange(minc, maxc+step, step)  # 1 cm steps
            
            # Four walls
            west  = np.column_stack([np.full_like(coords, minc), coords])
            east  = np.column_stack([np.full_like(coords, maxc), coords])
            south = np.column_stack([coords, np.full_like(coords, minc)])
            north = np.column_stack([coords, np.full_like(coords, maxc)])
            
            # Stack all together
            boundary_points[curenv] = np.vstack([west, east, south, north])           
      
            
        elif curenv == 'Hall':
           scaling = 4 #make length longer than width (vertical hallway)
           minx = -np.sqrt(area/scaling)/2
           maxx = np.sqrt(area/scaling)/2
           miny = minx*scaling
           maxy= maxx*scaling
           
           #Not a clever way to do this, but just need something that works for now
           westx = np.ones(np.ceil(maxy-miny).astype('int'))*minx
           westy = np.linspace(miny,maxy,np.ceil(maxy-miny).astype('int'))
           
           northx = np.linspace(minx,maxx,np.ceil(maxx-minx).astype('int'))
           northy = np.ones(np.ceil(maxx-minx).astype('int'))*maxy
           
           eastx = np.ones(np.ceil(maxy-miny).astype('int'))*maxx
           easty = np.linspace(maxy,miny,np.ceil(maxy-miny).astype('int'))
           
           southx = np.linspace(maxx,minx,np.ceil(maxx-minx).astype('int'))
           southy = np.ones(np.ceil(maxx-minx).astype('int'))*miny
           
           bounds_x = np.hstack([westx,northx,eastx,southx])
           bounds_y = np.hstack([westy,northy,easty,southy])
           
           boundpairs = np.vstack([bounds_x,bounds_y]).transpose()
           
           boundary_points[curenv] = boundpairs
            
        elif curenv == 'Circle':
             radius = np.sqrt(area/np.pi)
             circum = 2*np.pi*radius
             res = np.ceil(circum) #one boundary point for every 1 cm of wall
             thetas = np.linspace(0,2*np.pi,int(res))
             
             bounds_x = radius * np.cos(thetas)
             bounds_y = radius * np.sin(thetas)
             
             boundpairs = np.vstack([bounds_x,bounds_y]).transpose()
           
             boundary_points[curenv] = boundpairs
             
            
        else:
            
            print('Environment type not coded')

    return env_types, boundary_points

def set_trajstats():
    
    fs = 30 #sampling rate for video is often ~30 Hz. 
    
    'List of values from Raudies + Hasselmo 2012 for square and circular environment' 
    '0th values are for "Walking" '
    b_list = [6, 16.99,16.44,13.25,13.02] 
    mu_list = [0, -2.48,0.31,0.62,-0.03,1.89]
    sigma_list = [30, 350.58, 355.35,337.93,330.12,331.07]
    
    #inds = np.random.randint(0,3,(1,3)) #weirdly this generates an array of arrays so have to double index
    # b = b_list[inds[0][0]]
    # mu = mu_list[inds[0][1]]
    # sigma = sigma_list[inds[0][2]]

    b = b_list[0] 
    mu = mu_list[0]    #Used to be -2.48
    sigma = sigma_list[0] #Used to be 350.58 (?)

    return fs,b,mu,sigma

def analyzetraj(simx, simy, headingdir, boundary_points):
    #Collection of distances to all boundary points
    allodis = np.ones((np.size(simx), np.shape(boundary_points)[0]))
    
    #Distance to closest boundary point
    egodis = np.ones((np.size(simx),3)) # (distance,x,y)
    egoangle =np.ones(np.size(simx))

    #Angle to this closest boundary point
    posangl=np.ones(np.size(simx))
    closest_point =np.ones((np.size(simx), 2))
    
    
    for i in range(0,np.size(simx)): 
        
        curx = simx[i]
        cury = simy[i]
        
        allodis[i, :] = np.sqrt((curx - boundary_points[:,0])**2 + (cury - boundary_points[:,1])**2)
        closest_idx=np.argmin(allodis[i, :])
        closest_point[i] = boundary_points[closest_idx]
        egodis[i] = [allodis[i, closest_idx], closest_point[i,0], closest_point[i,1]]
        
        posangl[i] = np.arctan2(closest_point[i][1] - cury, closest_point[i][0] - curx) #angle from animal to closest boundary point.
        
        if posangl[i] < 0: #Always make this angle positive just for ease of debugging
            posangl[i] = np.mod(posangl[i], np.pi*2)
        
        egoangle[i] = np.mod(posangl[i]-headingdir[i],np.pi*2) #left is 3*pi/2 (not -pi/2)
      
    return allodis,egodis,egoangle,posangl, closest_point


"Running/Walking Policies"
def WallWalk(dur, fs, b, mu, sigma, boundary_points, curenv):
    "Simulates a rat walking back and forth along a single randomly chosen wall."
  
    dt = 1 / fs
    samples = int(dur * fs)
    
    simx = np.zeros(samples)
    simy = np.zeros(samples)
    headingdir = np.zeros(samples)
    
    speed = b #constant
    
    wall_choice = np.random.choice(['N', 'S', 'E', 'W'])
    
    min_x = np.min(boundary_points[curenv][:,0])
    max_x = np.max(boundary_points[curenv][:,0])
    min_y = np.min(boundary_points[curenv][:,1])
    max_y = np.max(boundary_points[curenv][:,1])
    
    # Initial positions, movement axis, and heading 
    if wall_choice == 'N': # Top wall
        y_const = max_y
        x_curr = min_x
        current_heading = 0 # Facing East initially
        moving_axis = 'X'
        direction_sign = 1
    elif wall_choice == 'S': # Bottom wall
        y_const = min_y
        x_curr = min_x
        current_heading = 0 # Facing East
        moving_axis = 'X'
        direction_sign = 1
    elif wall_choice == 'E': # Right wall
        x_const = max_x
        y_curr = min_y
        current_heading = np.pi/2 # Facing North initially
        moving_axis = 'Y'
        direction_sign = 1
    elif wall_choice == 'W': # Left wall
        x_const = min_x
        y_curr = min_y
        current_heading = np.pi/2 # Facing North
        moving_axis = 'Y'
        direction_sign = 1

    # Set the very first data point
    if moving_axis == 'X':
        simx[0] = x_curr
        simy[0] = y_const
    else:
        simx[0] = x_const
        simy[0] = y_curr
    headingdir[0] = current_heading
    
    # Run the back-and-forth simulation
    for i in range(1, samples):
        # Calculate how far the rat moves this frame
        delta = speed * dt * direction_sign
        
        if moving_axis == 'X':
            x_curr += delta
            
            # Check if rat hit the East corner
            if x_curr >= max_x:
                x_curr = max_x
                direction_sign = -1
                current_heading = np.pi # Turn around to face West
                
            # Check if rat hit the West corner
            elif x_curr <= min_x:
                x_curr = min_x
                direction_sign = 1
                current_heading = 0 # Turn around to face East
            
            simx[i] = x_curr
            simy[i] = y_const
            
        elif moving_axis == 'Y':
            y_curr += delta
            
            # Check if rat hit the North corner
            if y_curr >= max_y:
                y_curr = max_y
                direction_sign = -1
                current_heading = 3 * np.pi / 2 # Turn around to face South
                
            # Check if rat hit the South corner
            elif y_curr <= min_y:
                y_curr = min_y
                direction_sign = 1
                current_heading = np.pi / 2 # Turn around to face North
                
            simx[i] = x_const
            simy[i] = y_curr
        
        headingdir[i] = current_heading
        
    return simx, simy, headingdir

def make_envirs(numenv): 
    
    possible_types = ['Square','Hall','Circle',]
    the_types = possible_types[0:numenv]
    #For now this is hard coded to just do square, then hall, then circle
    #only writing the square and hall code for now for ease
    env_types = []
    boundary_points = dict()
    
    area = 10000 #assuming environment as area of 10,000 cm^2
    #so for square environment have a 1 m x 1 m box
    # setting center of the environment to be approximately 0,0
    #Using ~1 cm resolution for boundary points
    
    for curenv in the_types:
        print(curenv)
        env_types.append(curenv)
        
        if curenv == 'Square': #basically produce an array that contains all the points along the walls
            #Use 1 cm bins for each boundary point
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
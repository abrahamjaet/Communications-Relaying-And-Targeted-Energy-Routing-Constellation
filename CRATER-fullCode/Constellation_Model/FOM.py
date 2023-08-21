import matplotlib.pyplot as plt
import numpy as np


moon_rad = 1740 # radius of the Moon in [km]
lat_spacing = 10 # [deg]
max_points = 36 
latitude = np.linspace(90,-90,19)
longitude = np.linspace(0,350,36)


# print(longitude)

# print(latitude)

points = []
counter = 0
# northern hemisphere 
for i in latitude: 
    radians = np.radians(i)
    points += [abs(round(max_points*np.cos(radians)))]
    counter = counter+1
print(points)

ax = plt.axes(projection = '3d')
coordinates = []
i = 0
for point in points: 
    longitude = np.linspace(0,350,point)
    lat = latitude[i]
    
    for long in longitude:
        x = moon_rad*np.cos(np.radians(lat))*np.cos(np.radians(long))
        y = moon_rad*np.cos(np.radians(lat))*np.sin(np.radians(long))
        z = moon_rad*np.sin(np.radians(lat))
        data = [x,y,z,lat,long]
        ax.scatter3D(x,y,z)
        coordinates.append(data)
    i = i + 1 
# print(coordinates)
test = coordinates[0]
test = test[0]
print(test)
print(coordinates[0])
print(coordinates[1])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

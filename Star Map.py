'''
Code to create a sky map of the stars
    in the Healy-Horiuchi 23 catalog
'''

import matplotlib
from matplotlib import pyplot as plt
import scienceplots
from IPython import display
plt.style.use(["science","ieee","bright"])
import numpy as np
from numpy import pi as pi
from numpy import sqrt, cos, sin, exp, log
from matplotlib import ticker, cm
from matplotlib import animation
from matplotlib.animation import FuncAnimation
import math
import scipy as sp
from scipy import special as spc


months_to_eq_days = {"January":287, "February":318,
                     "March":-19,"April":12,
                     "May":42, "June":73,
                     "July":103, "August":134,
                     "September":165,"October":195,
                     "November":226, "December":256}

city_to_latitude = {"New York":40.75, "New Orleans":30,
                    "Blacksburg":37, "Los Angeles":34,
                    "Chicago":41.8,"Atlanta":33.75,
                    "Charolette":35.25, "Houston":29.75,
                    "Seatle":47.6, "Louisville":38.25,
                    "Nashville": 36.16}

while True:
    month = input("Input a Month: ")
    try:
        months_to_eq_days[month]
        break
    except:
        print("Please Enter Valid Month")
        print(months_to_eq_days.keys())

while True:
    day = input("Input day of month: ")
    try:
        int(day)
        break
    except:
        print("Please Enter Valid Day")
        
while True:
    time = input("Input time (in hours) after midnight: ")
    try:
        if float(time) >0 and float(time) < 24:
            break
        else:
            print("Please Enter Valid Time")
    except:
        print("Please Enter Valid Time")

while True:
    city = input("Input a City: ")
    try:
        city_to_latitude[city]
        break
    except:
        print("Please Enter A Given City")
        print(city_to_latitude)




days_after_equinox = months_to_eq_days[month] + int(day)
rev_angle = days_after_equinox * (2*pi/365)

hours_after_midnight = float(time)
rot_angle = hours_after_midnight * (pi/12)

tot_rot_angle = (-rev_angle - rot_angle) + pi

latitude_deg = city_to_latitude[city]
lat_rad = latitude_deg * (pi/180)

tilt_rad = 22.5 * (pi/180)

starting_direction = np.array([cos(lat_rad),0,sin(lat_rad)])
starting_north = np.array([0,0,1])

rev_rot = np.array([[cos(-tot_rot_angle), -sin(-tot_rot_angle),0],
                    [sin(-tot_rot_angle), cos(-tot_rot_angle),0],
                    [0,0,1]])

direction = np.matmul(rev_rot,starting_direction)

int_north = np.matmul(rev_rot,starting_north)

tilt_rot = np.array([[1,0,0],
                     [0, cos(tilt_rad), sin(tilt_rad)],
                     [0, -sin(tilt_rad),cos(tilt_rad)]])

ec_direction = np.matmul(tilt_rot,direction)
end_north = np.matmul(tilt_rot,int_north)

north_dir = end_north - np.sum(end_north*ec_direction)*ec_direction

north_dir = north_dir / sqrt(sum(north_dir*north_dir))

east_dir = -np.array([ec_direction[1]*north_dir[2] - ec_direction[2]*north_dir[1],
                     ec_direction[2]*north_dir[0] - ec_direction[0]*north_dir[2],
                     ec_direction[0]*north_dir[1] - ec_direction[1]*north_dir[0]])
#Example Zodiac Constelations



Constellations_angs = {"Ares":[28,53],
                           "Taurus":[53,90],
                           "Gemini":[90,118],
                           "Cancer":[118,138],
                           "Leo":[138,174],
                           "Virgo":[174,218],
                           "Libra":[218,241],
                           "Scorpius":[241,248],
                           "Ophiuchun":[248,266],
                           "Sagittarius":[266,300],
                           "Capricornus":[300,327.5],
                           "Aquarius":[327.5,351.7],
                           "Pisces":[351.7,360 + 28.687]}

fig = plt.figure()
x_vals = np.linspace(-1,1,1000)
plt.fill_between(x_vals,sqrt(1-x_vals**2),-sqrt(1-x_vals**2),color = "Black",alpha = 0.5)

Opening_angle = pi/2
open_cos = cos(Opening_angle)

for key in Constellations_angs.keys():
    angs = Constellations_angs[key]
    
    pos = np.array([cos(np.mean(angs)*pi/180),sin(np.mean(angs)*pi/180),0])
    
    #print(key)
    dot_prod = sum(pos*ec_direction)
    #print(dot_prod)
    
    if dot_prod > open_cos:
        print(key)
        #print("north",sum(north_dir*pos))
        #print("east", sum(east_dir*pos))
        
        phi = np.arccos(sum(-east_dir*pos)) * np.sign(sum(north_dir*pos))
        #phi = np.arcsin(sum(north_dir*pos)) + pi * np.heaviside(sum(-east_dir*pos),0)
        sin_theta = sqrt(1-dot_prod**2)
        
        plt.scatter(sin_theta*cos(phi), sin_theta*sin(phi), marker = "*", label = key)
        
#plt.legend()
plt.xlim([-1,1])
plt.ylim([-1,1])
    


star_pos = np.zeros((674,3))
file = open("Star_Catalog_Lon_Lat.csv")
line_index = 0
for line in file:
    line = line.split(",")
    if line_index < 1:
        #print(line[0])
        nothing = 0
    else:
        lon = float(line[0])*pi/180 
        lat = float(line[1])*pi/180
        ind_pos = np.array([cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)])
        star_pos[line_index-1] = ind_pos
    line_index += 1

file.close()
    
 
fig = plt.figure()
x_vals = np.linspace(-1,1,1000)
plt.fill_between(x_vals,sqrt(1-x_vals**2),-sqrt(1-x_vals**2),color = "Black",alpha = 1)

Opening_angle = pi/2
open_cos = cos(Opening_angle)    

for pos in star_pos:
    dot_prod = sum(pos*ec_direction)
    if dot_prod > open_cos:
        #print("north",sum(north_dir*pos))
        #print("east", sum(east_dir*pos))
        
        phi = np.arccos(sum(-east_dir*pos)) * np.sign(sum(north_dir*pos))
        #phi = np.arcsin(sum(north_dir*pos)) + pi * np.heaviside(sum(-east_dir*pos),0)
        sin_theta = sqrt(1-dot_prod**2)
        
        plt.scatter(sin_theta*cos(phi), sin_theta*sin(phi),s=1, marker = "*",color = "red")
        
fig2 = plt.figure()

x_vals = np.linspace(-1,1,1000)
plt.fill_between(x_vals,sqrt(1-x_vals**2),-sqrt(1-x_vals**2),color = "Black")

points_multi_plotted = plt.plot([],color = "Black",marker = "*", mfc = "red", mec = "red",markersize = 1)

points_plotted = points_multi_plotted[0]

figtext = plt.text(.8,.9,city)
def AnimationFunction(frame):
    hours_after_midnight = 18 + 24*frame/100
    
    mil_hour = round((hours_after_midnight % 24)//1)
    mil_minute = round((hours_after_midnight % 1) * 60)
    rot_angle = hours_after_midnight * (pi/12)
    

    tot_rot_angle = (-rev_angle - rot_angle) + pi

    latitude_deg = city_to_latitude[city]
    lat_rad = latitude_deg * (pi/180)

    tilt_rad = 22.5 * (pi/180)

    starting_direction = np.array([cos(lat_rad),0,sin(lat_rad)])
    starting_north = np.array([0,0,1])

    rev_rot = np.array([[cos(-tot_rot_angle), -sin(-tot_rot_angle),0],
                        [sin(-tot_rot_angle), cos(-tot_rot_angle),0],
                        [0,0,1]])

    direction = np.matmul(rev_rot,starting_direction)

    int_north = np.matmul(rev_rot,starting_north)

    tilt_rot = np.array([[1,0,0],
                         [0, cos(tilt_rad), sin(tilt_rad)],
                         [0, -sin(tilt_rad),cos(tilt_rad)]])

    ec_direction = np.matmul(tilt_rot,direction)
    
    end_north = np.matmul(tilt_rot,int_north)

    north_dir = end_north - np.sum(end_north*ec_direction)*ec_direction

    north_dir = north_dir / sqrt(sum(north_dir*north_dir))

    east_dir = -np.array([ec_direction[1]*north_dir[2] - ec_direction[2]*north_dir[1],
                         ec_direction[2]*north_dir[0] - ec_direction[0]*north_dir[2],
                         ec_direction[0]*north_dir[1] - ec_direction[1]*north_dir[0]])
    
    xstars = np.array([])
    ystars = np.array([])
    
    Opening_angle = pi/2
    open_cos = cos(Opening_angle)    

    for pos in star_pos:
        dot_prod = sum(pos*ec_direction)
        if dot_prod > open_cos:
            #print("north",sum(north_dir*pos))
            #print("east", sum(east_dir*pos))
            
            phi = np.arccos(sum(-east_dir*pos)) * np.sign(sum(north_dir*pos))
            #phi = np.arcsin(sum(north_dir*pos)) + pi * np.heaviside(sum(-east_dir*pos),0)
            sin_theta = sqrt(1-dot_prod**2)
            
            xstars = np.append(xstars,sin_theta*cos(phi))
            ystars = np.append(ystars, sin_theta*sin(phi))
        
    points_plotted.set_data(xstars,ystars)
    #figtext.set_text("t = "+str(mil_hour) +":" +str(mil_minute))
    
anim_created = FuncAnimation(fig2, AnimationFunction,frames = 100,interval = 25)

# To save the animation using Pillow as a gif
writer = animation.PillowWriter(fps=15,
                                 metadata=dict(artist='Me'),
                                 bitrate=1800)
anim_created.save(city+'.gif', writer=writer)
'''
video = anim_created.to_html5_video()
html = display.HTML(video)
display.display(html)
'''
'''
fig3 = plt.figure()

x_vals = np.linspace(-1,1,1000)
plt.fill_between(x_vals,sqrt(1-x_vals**2),-sqrt(1-x_vals**2),color = "Black")

points_multi_plotted_2 = plt.plot([],color = "Black",marker = "*", mfc = "red", mec = "red",markersize = 1)

points_plotted_2 = points_multi_plotted_2[0]

figtext = plt.text(.8,.9,"Annual Change")
def AnimationFunction2(frame):
    
    days_after_equinox = frame * 365/100
    rev_angle = days_after_equinox * (2*pi/365)
    
    hours_after_midnight = float(time)
    
    rot_angle = hours_after_midnight * (pi/12)
    
    tot_rot_angle = (-rev_angle - rot_angle) + pi

    latitude_deg = city_to_latitude[city]
    lat_rad = latitude_deg * (pi/180)

    tilt_rad = 22.5 * (pi/180)

    starting_direction = np.array([cos(lat_rad),0,sin(lat_rad)])
    starting_north = np.array([0,0,1])

    rev_rot = np.array([[cos(-tot_rot_angle), -sin(-tot_rot_angle),0],
                        [sin(-tot_rot_angle), cos(-tot_rot_angle),0],
                        [0,0,1]])

    direction = np.matmul(rev_rot,starting_direction)

    int_north = np.matmul(rev_rot,starting_north)

    tilt_rot = np.array([[1,0,0],
                         [0, cos(tilt_rad), sin(tilt_rad)],
                         [0, -sin(tilt_rad),cos(tilt_rad)]])

    ec_direction = np.matmul(tilt_rot,direction)
    end_north = np.matmul(tilt_rot,int_north)

    north_dir = end_north - np.sum(end_north*ec_direction)*ec_direction

    north_dir = north_dir / sqrt(sum(north_dir*north_dir))

    east_dir = -np.array([ec_direction[1]*north_dir[2] - ec_direction[2]*north_dir[1],
                         ec_direction[2]*north_dir[0] - ec_direction[0]*north_dir[2],
                         ec_direction[0]*north_dir[1] - ec_direction[1]*north_dir[0]])
    
    xstars2 = np.array([])
    ystars2 = np.array([])
    
    Opening_angle = pi/2
    open_cos = cos(Opening_angle)    
    
    included = 0
    for pos in star_pos:
        dot_prod = sum(pos*ec_direction)
        if dot_prod > open_cos:
            included += 1
            #print("north",sum(north_dir*pos))
            #print("east", sum(east_dir*pos))
            
            phi = np.arccos(sum(-east_dir*pos)) * np.sign(sum(north_dir*pos))
            #phi = np.arcsin(sum(north_dir*pos)) + pi * np.heaviside(sum(-east_dir*pos),0)
            sin_theta = sqrt(1-dot_prod**2)
            
            xstars2 = np.append(xstars2,sin_theta*cos(phi))
            ystars2 = np.append(ystars2, sin_theta*sin(phi))
        
    points_plotted_2.set_data(xstars2,ystars2)
    print('stars',included)
    #figtext.set_text("t = "+str(mil_hour) +":" +str(mil_minute))
    
anim_created2 = FuncAnimation(fig3, AnimationFunction2,frames = 100,interval = 25)

video2 = anim_created2.to_html5_video()
html2 = display.HTML(video2)
display.display(html2)
    
'''

# Information Technology in Ship Design and Production
"""
If this program is running for the first time in the system which doesn't
have the reqired libraries already installed in it. It will automaticallly
install them from the online python library. It is advisable for the system to
be connected to the internet.

This program may not work properly in python ver. 2.7. It is advisable to use
spyder (cross-platform for Python) which will have majority of the libraries
installed in it.

Please make sure the "STL" file and the python file is present in the same folder
before executing the program. The output will be saved as 'output.csv' file and
diagrams are saved in pdf format in the same folder.

The user has to give the input in the 'config.ini' which also be placed in the
same folder as the code. Default input has been given in the config file, furthur
changes can be made by the user. Specify the unit of the stl file in the config file
since stl files does not save the unit informations. 'Meter' is recommended.

The Input STL file must contain only one symmetric half of the ship Hull, but the output
will be calculated for complete Hull. Input Hull with positive coordinates would be
advisable for easy visualization. Please ignore any error message regarding the STL file
format.
"""
import pip
import numpy as np
from stl import mesh
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import trimesh
import math
import configparser
import csv

import concurrent.futures

config = configparser.ConfigParser()
config.read("config.ini")
config.sections()

def ConfigSectionMap(section):
    dict1 = {}
    options = config.options(section)
    for option in options:
        try:
            dict1[option] = config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

#Entering the name of the stl file
fname = ConfigSectionMap("STLtoHydrostatics")['file_name']
m = mesh.Mesh.from_file(fname)
#np.savetxt("mesh.csv", m, delimiter=",")
#print(m[0])


unit = str(ConfigSectionMap("STLtoHydrostatics")['units(m_or_mm)'])

# The mesh normals
m.normals
# Declaring the mesh vectors
m.v0, m.v1, m.v2
# Storing the vertices(v1, v2, v3) into 3 mesh vectors(v0, v1, v2).
assert (m.points[0][0:3] == m.v0[0]).all()
assert (m.points[0][3:6] == m.v1[0]).all()
assert (m.points[0][6:9] == m.v2[0]).all()
assert (m.points[1][0:3] == m.v0[1]).all()


#plt.xlim(5000, 12500)
#plt.ylim(-3750, 3750)
#plt.zlim(-100, 1000)

#Merging all the vectors into one
v = np.concatenate((m.v0, m.v1, m.v2), axis=0)
np.savetxt("num0.csv", m.v0, delimiter=",")
np.savetxt("num1.csv", m.v1, delimiter=",")
np.savetxt("num2.csv", m.v2, delimiter=",")
np.savetxt("num.csv", v, delimiter=",")
num_ver = v.size/3
num_tri = num_ver/3
print("number of vertices =", num_ver)
print("number of triangles =", num_tri)

tri_array = np.empty((0, 3))
with concurrent.futures.ProcessPoolExecutor() as executor:
    for triangle in range(0, int(num_tri), 1):
        tri = m[triangle].reshape((3,3))
        tri_array = np.append(tri_array,tri, axis=0)


#np.savetxt("tri.csv", tri_array, delimiter=",")
#print(tri_array)
np.set_printoptions(threshold=np.nan)


#Class is created to calculate length, height and breadth of the bounding box/Hull   

def length(v):
    x= v[:,0]
    a = x.min()
    b = x.max()
    l = b-a
    return(l)
    
def breadth(v):
    y= v[:,1]
    a = y.min()
    b = y.max()
    l = 2*(b-a)
    return(l)
    
def height(v):
    z= v[:,2]
    a = z.min()
    b = z.max()
    l = b-a
    return(l)
    
print('\nThe overall Length(x axis) of the hull in %s is: '%unit)
print(length(v))
print('The overall Breadth(y axis) of the hull in %s is: '%unit)
print(breadth(v))
print('The overall Height(z axis) of the hull in %s is: '%unit)
print(height(v))

def sortpts(inter):
    inter_yp = np.empty((0, 3))
    inter_ym = np.empty((0, 3))
    inter_y1 = np.empty((0, 3))
    inter_y2 = np.empty((0, 3))
    inter_y0 = np.empty((0, 3))
    lim = int(inter.size/3)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0, lim, 1):
            a = inter[:,1][i]
            if a > 0:
                y1 = np.array([inter[i]])
                inter_y1 = np.append(inter_y1, y1, axis = 0)
            elif a == 0:
                y0 = np.array([inter[i]])
                inter_y0 = np.append(inter_y0, y0, axis = 0)
            elif a < 0:
                y2 = np.array([inter[i]])
                inter_y2 = np.append(inter_y2, y2, axis = 0)
    inter_yp = inter_y1[inter_y1[:,0].argsort()]
    inter_yn = inter_y0[inter_y0[:,0].argsort()]
    inter_ym = inter_y2[(-inter_y2[:,0]).argsort()]
    inpts = np.concatenate((inter_yp, inter_yn, inter_ym), axis = 0)
    return(inpts)

def point_plane(z0, z1, plane):
    d0 = plane - z0
    d1 = plane - z1
    if d0*d1 < 0:
        p = d0/(d0-d1)
        #print(t)
        return(p)
    elif d0*d1 == 0:
        p = 0
        return(p)
    else:
        p = str('false')
        return(p)

def point_point(x0, y0, x1, y1):
    d = math.sqrt(math.pow(x1 - x0, 2) + math.pow(y1 - y0, 2))
    return(d)

def moment_area(inter_pt):
    lim = int((inter_pt.size / 3)-1)
    area_trap =0
    i_x = 0
    cent_trap = 0
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0, lim, 1):
            h = abs(inter_pt[:,0][i] - inter_pt[:,0][i+1])
            a = abs(inter_pt[:,1][i])*2# - inter_pt[:,1].min())*2
            b = abs(inter_pt[:,1][i+1])*2# - inter_pt[:,1].min())*2

            area_trap += ((a+b)/2)*h
            i_x += (h/48)*(a+b)*(math.pow(a, 2) + math.pow(b, 2))
            cent_trap += (((a+b)/2)*h)*(abs(inter_pt[:,0][i])+(h/2))
            continue    
        
    return(area_trap, i_x, cent_trap)

def volume_wp(array):
    lim =int(((array.size)/2)-1)
    vol_trap = 0
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for i in range(0, lim, 1):
            h = abs(array[:,0][i+1] -array[:,0][i])
            a = abs(array[:,1][i])
            b = abs(array[:,1][i+1])
            c = abs(b - a)
            
            vol_trap += (a*h)+(0.5*h*c)
            continue
        
    return(vol_trap)

def vol_immersed(inter_pt, tri_array, d):
    lim = int(tri_array.size/3)
    vol = np.empty((0, 3))
    for i in range(0, lim, 1):
        b = tri_array[:, 2][i]
        if b < d:
            all_ver = np.array([tri_array[i]])
            vol = np.append(vol, all_ver, axis =0)
    ver = np.concatenate((vol, inter_pt), axis =0)
    return(ver)

def block_coef(volume, inter_pt, draft):
    l = inter_pt[:,0].max() - inter_pt[:,0].min()
    b = inter_pt[:,1].max() - inter_pt[:,1].min()
    if draft == 0:
        c_b=0
    else:
        c_b = volume/(l*b*draft)
    return(c_b)
    
def center_grav(moment, volume, gm):
    if volume == 0:
        kg = 0
    else:
        kg = (moment/volume) - gm
        #print(kg)
    return(kg)


    
def main():
    d_e = float(ConfigSectionMap("STLtoHydrostatics")['max_draft'])
    sec = int(ConfigSectionMap("STLtoHydrostatics")['no_of_waterplanes'])
    d_st = float(d_e/sec)
        
    gm = float(ConfigSectionMap("STLtoHydrostatics")['gm_dist'])

    centroid_arr = np.empty((0, 1))
    area_wp = np.empty((0, 1))
    area_arr= np.empty((0, 2))
    vol_arr = np.empty((0, 1))
    mom_arr = np.empty((0, 1))
    cb_arr = np.empty((0, 1))
    cog = np.empty((0,1))
    draft_arr = np.empty((0,1))
    len_arr = np.empty((0, 1))
    bre_arr = np.empty((0, 1))
    
    x_buo =0
    z_buo =0

    buo_arr = np.empty((0, 3))


    with concurrent.futures.ProcessPoolExecutor() as executor:    
        for draft in np.arange(0, d_e + d_st, d_st):
            z= v[:,2]
            d = z.min() + draft

            
            draft_d = np.array([draft]).reshape(1,1)
            draft_arr = np.append(draft_arr, draft_d, axis=0)
            
            inter_l1 = np.empty((0, 3))
            inter_l2 = np.empty((0, 3))
            inter_l3 = np.empty((0, 3))
            inter_pt = np.empty((0, 3))

            with concurrent.futures.ProcessPoolExecutor() as executor:
                for i in range(0, int(num_ver), 3):
                            
                    t = point_plane(tri_array[:,2][i], tri_array[:,2][i+1], d)
                    if t!='false':
                        r_t = np.array([tri_array[i] + (t* (tri_array[i+1] - tri_array[i]))])
                        inter_l1 = np.append(inter_l1, r_t, axis=0)

                    t = point_plane(tri_array[:,2][i+1], tri_array[:,2][i+2], d)
                    if t!='false':
                        r_t = np.array([tri_array[i+1] + (t* (tri_array[i+2] - tri_array[i+1]))])
                        inter_l2 = np.append(inter_l2, r_t, axis=0)

                    t = point_plane(tri_array[:,2][i+2], tri_array[:,2][i], d)
                    if t!='false':
                        r_t = np.array([tri_array[i+2] + (t* (tri_array[i] - tri_array[i+2]))])
                        inter_l3 = np.append(inter_l3, r_t, axis=0)
                    continue
                
            inter = np.concatenate((inter_l1, inter_l2, inter_l3), axis=0)
            inter_pt = sortpts(inter)

            v_imer = vol_immersed(inter_pt, tri_array, d)
            l = length(v_imer)
            b = breadth(v_imer)

            len_arr = np.append(len_arr, np.array([l]).reshape(1,1), axis = 0)
            bre_arr = np.append(bre_arr, np.array([b]).reshape(1,1), axis = 0)

            area, moment, cent_wp = moment_area(inter_pt)

            area_wp = np.append(area_wp, np.array([area]).reshape(1,1), axis =0)
            arr = np.array([d, area]).reshape(1, 2)
            area_arr = np.append(area_arr, arr, axis =0)
            
            cent_wp = cent_wp/area
            centroid = np.array([cent_wp]).reshape(1,1)
            centroid_arr = np.append(centroid_arr, centroid, axis=0)

            
            mom_wp = np.array([moment]).reshape(1, 1)
            mom_arr = np.append(mom_arr, mom_wp, axis=0)
            
            volume = volume_wp(area_arr)
                
            vol = np.array([volume]).reshape(1,1)
            vol_arr = np.append(vol_arr, vol, axis=0)

            a_wp = area
            x_wp = cent_wp
            z_wp = draft
                                    
            x_buo += a_wp*x_wp        
            z_buo += a_wp*z_wp

            buo = np.array([x_buo, 0, z_buo]).reshape(1, 3)
            buo_arr = np.append(buo_arr, buo, axis=0)

            if draft == d_e:
                lim = int((buo_arr.size/3))
                center_of_buo_x = np.empty((0, 1))
                center_of_buo_z = np.empty((0, 1))
                with concurrent.futures.ProcessPoolExecutor() as executor:
                   for i in range (0, lim, 1):
                        area_sum = np.cumsum(area_arr[:,1])
                        center_of_buo_x = np.append(center_of_buo_x, np.array([buo_arr[:,0][i]/area_sum[i]]).reshape(1, 1), axis=0)
                        center_of_buo_z = np.append(center_of_buo_z, np.array([buo_arr[:,2][i]/area_sum[i]]).reshape(1, 1), axis=0)
                        continue

            cb = block_coef(volume, inter_pt, draft)
            c_b = np.array([cb]).reshape(1,1)
            cb_arr = np.append(cb_arr, c_b, axis= 0)
            
            g = center_grav(moment, volume, gm)
            k_g = np.array([g]).reshape(1,1)
            cog = np.append(cog, k_g, axis=0)
            continue
    cog_arr = np.empty((0,1))
    for i in range (0, int((cog.size)), 1):
        cog_arr = np.append(cog_arr, np.array(cog[:,0][i]+center_of_buo_z[:,0][i]).reshape(1,1), axis=0)
        continue

    output = np.concatenate((draft_arr, vol_arr, cb_arr, center_of_buo_z, area_wp, mom_arr, cog_arr, len_arr, bre_arr), axis=1)
    np.savetxt("Output.csv", output, delimiter=",")
    
    np.savetxt("center of gravity.csv", cog_arr, delimiter=",")
    np.savetxt("Block Coeffecient.csv", cb_arr, delimiter=",")    
    np.savetxt("Center of Buoyance x.csv", center_of_buo_x, delimiter=",")
    np.savetxt("Center of Buoyance z.csv", center_of_buo_z, delimiter=",")
    np.savetxt("Moment of Inertia wp.csv", mom_arr, delimiter=",")
    np.savetxt("Centroid.csv", centroid_arr, delimiter=",")
    np.savetxt("Area wp.csv", area_wp, delimiter=",")
    np.savetxt("Area.csv", area_arr, delimiter=",")
    np.savetxt("Volume.csv", vol_arr, delimiter=",")  
    print(volume)


    with open('output.csv',newline='') as f:
        r = csv.reader(f)
        data = [line for line in r]
    with open('output.csv','w',newline='') as f:
        w = csv.writer(f)
        w.writerow(['Draft','Displacement','C_b', 'KB', 'A_w', 'I_y', 'KG', 'Length', 'Breadth'])
        w.writerows(data)

    fig1, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Height in %s'%unit)
    plt.title('Height of Center of Buoyancy and Gravity from base of ship')
    ax.plot(draft_arr, center_of_buo_z, 'r--', label='Height of COB in %s'%unit)
    ax.plot(draft_arr, cog_arr, 'b:', label='Height of COG in %s'%unit)

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig1.savefig('KB,KG.pdf')

    fig2, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Length, Breadth in %s'%unit)
    plt.title('Length and Breadth of Displaced volume for each draft')
    ax.plot(draft_arr, len_arr, 'r--', label='Length  in %s'%unit)
    ax.plot(draft_arr, cog_arr, 'b:', label='Breadth in %s'%unit)

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig2.savefig('len,bre.pdf')

    fig3, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Volume in $%s^3$'%unit)
    plt.title('Displacement curve based on the draft range')
    ax.plot(draft_arr, vol_arr, 'b', label='Displacement in $%s^3$'%unit)

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig3.savefig('Disp.pdf')

    fig4, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Area in $%s^2$'%unit)
    plt.title('Waterplane area for each draft')
    ax.plot(draft_arr, area_wp, 'b', label='Waterplane area in $%s^2$'%unit)

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig4.savefig('area_wp.pdf')

    fig5, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Moment of Inertia $%s^4$'%unit)
    plt.title('Moment of Inertia about x-axis for each waterplane')
    ax.plot(draft_arr, mom_arr, 'b', label='Moment of Inertia about x-axis $%s^4$'%unit)

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig5.savefig('moment of inertia.pdf')
    
    fig6, ax = plt.subplots()
    plt.xlabel('Draft in %s'%unit)
    plt.ylabel('Block Coefficient (Dimensionless)')
    plt.title('Block Coefficient $C_B$ for each draft')
    ax.plot(draft_arr, mom_arr, 'b', label='$C_B$')

    legend = ax.legend(loc='upper left', shadow=False)
    for label in legend.get_texts():
        label.set_fontsize('medium')

    for label in legend.get_lines():
        label.set_linewidth(1.5)
    plt.show()
    fig6.savefig('Block Coefficient.pdf')

if __name__ == '__main__':
    main()

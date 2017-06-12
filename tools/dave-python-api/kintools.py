def triangleList(faces, vertices):

    # init output string holder
    output = []
    
    # kinemage headers 
    output.append("@trianglelist {triangles}")

    for f in faces:

        # create a line for this triangular face
        line = "{}X"
        [v1, v2, v3] = f[0:3]
        for vert in [v1, v2, v3]:
            x,y,z = vertices[vert]
            line += " %f %f %f" %(x,y,z)
        
        # append the line to the output
        output.append(line)
        
    return "\n".join(output)

#
# Kinemage colour definitions
#
saturated_colours = ["red", "orange", "gold", "yellow", "lime", 
                    "green", "sea", "cyan", "sky", "blue",
                    "purple", "magenta", "hotpink"]
semi_saturated_colours = ["pink", "peach", "yellow", "sea", "sky", "lilac"]
pastel_colours = ["pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"]
neutral_colours = ["white", "gray", "brown"]

# make a list of all unique colours
all_colours = []
for colourlist in [saturated_colours, semi_saturated_colours, pastel_colours, neutral_colours]:
    for colour in colourlist:
        if colour not in all_colours:
            all_colours.append(colour)

# don't use these other colours
#other_colours = ["deadwhite", "deadblack", "invisible"]

def randomColour():
    from random import choice    
    return choice(all_colours)

def red_blue_colours():
    from math import floor, ceil
    
    output = []
    blue_hue = 240
    red_hue = 0
    underflow_hue = 80
    overflow_hue = 160
    sat = 100
    
    for ii in range(0,101):
        output.append("@hsvcolour {red_%d} %d %d %d" %(ii, red_hue, ii, sat))
        output.append("@hsvcolour {blue_%d} %d %d %d" %(ii, blue_hue, ii, sat))
    output.append("@hsvcolour {red_101} %d %d %d" %(underflow_hue, 100, sat))
    output.append("@hsvcolour {blue_101} %d %d %d" %(overflow_hue, 100, sat))
    
    return "\n".join(output)

def colour_ranges(num_ranges):
    from math import floor, ceil
    
    output = []
    colour_band_width = 360.0 / num_ranges
    colour_increment = 5.0
    num_colours_per_band = int(floor(colour_band_width / colour_increment))
    #ctr_list = []
    colour_names = []
    
    for ii in range(0,num_ranges):
        colour_names.append([])
        ctr = 0
        hstart = int(floor(ii*colour_band_width))
        hend = max(360,int(ceil((ii+1)*colour_band_width)))
        for inc_ctr in range(0, num_colours_per_band):
            hue = hstart + floor(inc_ctr * colour_increment)
            for saturation in range(35,101,5):
                colour_name = """Range%dColour%d""" %(ii, ctr)
                colour_names[ii].append(colour_name)
                output.append("@hsvcolor {%s} %d %d 100" %(colour_name, hue, saturation))
                ctr += 1
        #ctr_list.append(ctr)
    
    return colour_names, "\n".join(output)

def randomPairPastels():

    #assert(len(semi_saturated_colours) == len(pastel_colours))
    max_colour_idx = min([len(semi_saturated_colours), len(pastel_colours)])
    
    semi_saturated_colours
    pastel_colours
    from random import randint
    i = randint(0, max_colour_idx)
    
    return semi_saturated_colours[i], pastel_colours[i]

    
    
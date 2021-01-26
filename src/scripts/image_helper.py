import sys
import glob
import logging
import math
import time
import networkx as nx
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont
import shutil

sys.path.append('../../')
from config import *
sys.path.append(scripts_dir)
from utils import *

def cleanup_images(superimposition_output_dir, progressive_dir, representative_dir, subfamilies_dir):
    image_filelist = glob.glob(os.path.join(superimposition_output_dir, '*.png'))
    create_directory(progressive_dir)
    for file in image_filelist:
        out = os.path.join(progressive_dir, os.path.basename(file))
        img = Image.open(file)
        img = crop_image(img)
        img.save(out)
        img.close()
        os.remove(file)

    image_filelist = glob.glob(os.path.join(representative_dir, '*.png'))
    for file in image_filelist:
        img = Image.open(file)
        img = crop_image(img)
        img.save(file)
        img.close()

    create_directory(subfamilies_dir)
    combined_image_dir = os.path.join(superimposition_output_dir, 'subfamily', 'Subfamily_Combined')
    image_filelist = glob.glob(os.path.join(combined_image_dir, '*.png'))
    for file in image_filelist:
        shutil.copyfile(file, os.path.join(subfamilies_dir, os.path.basename(file)))
        # os.remove(file)
    # remove_all_from_dir(os.path.join(superimposition_output_dir, 'subfamily'))

def separate_family_name_and_loop_id(load_id):
    pcs = load_id.split('_')
    family_name = '_'.join(pcs[:-1])
    loop_id = pcs[-1]
    return family_name, loop_id

def get_node_name(node, allowed_characters = 6):
    # name, id = node.strip().split("_")
    name, id = separate_family_name_and_loop_id(node)
    if len(node) > allowed_characters:
        if len(name) > allowed_characters:
            if len(node) > allowed_characters * 2:
                dots = "..."
                id_len = len(id)
                name_len = allowed_characters * 2 - id_len - 4
                upper_dot_len = allowed_characters - name_len
                new_name = name[:name_len] + dots[:upper_dot_len] + "\n" + dots[upper_dot_len:] + "_" + id
            else:
                new_name = node[:allowed_characters] + "\n" + node[allowed_characters:]
        else:
            if len(name) < allowed_characters:
                new_name = name + "_\n" + id
            else:
                new_name = name + "\n_" + id
    else:
        new_name = node

    return new_name

def get_border_region(fixed_positions):
    x_min = y_min = x_max = y_max = 0.0
    for node in fixed_positions:
        x, y = fixed_positions[node]
        x_min = min(x_min, x)
        y_min = min(y_min, y)
        x_max = max(x_max, x)
        y_max = max(y_max, y)

    return x_min, y_min, x_max, y_max

def get_coordinates_in_circle(n):
    return_list = []
    for i in range(n):
        theta = float(i)/n*2*math.pi
        x = np.cos(theta)
        y = np.sin(theta)
        return_list.append((x,y))
    return return_list

def get_n_fixed_points_next_to_circle(n, outside_angle, point1, point2, scale=0.8):
    x1, y1 = point1
    x2, y2 = point2
    return_list = []
    for i in range(n):
        theta = -1 * (i+1) * outside_angle / (n+1)
        x = (x2-x1) * np.cos(theta) - (y2-y1) * np.sin(theta)
        y = (x2-x1) * np.sin(theta) + (y2-y1) * np.cos(theta)
        new_x = x+x1
        new_y = y+y1
        new_x *= scale
        new_y *= scale
        return_list.append((new_x,new_y))
    return return_list

def get_n_fixed_points_far_from_circle(n, point1, point2, scale=1):
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = (-1 * (y2-y1)+x1, (x2-x1)+y1)  # counter-clockwise perpendicular

    outside_angle = math.pi
    return_list = []
    for i in range(n):
        theta = (i+1) * outside_angle / (n+1)
        x = (x3-x1) * np.cos(theta) - (y3-y1) * np.sin(theta)
        y = (x3-x1) * np.sin(theta) + (y3-y1) * np.cos(theta)
        new_x = x+x1
        new_y = y+y1
        new_x *= scale
        new_y *= scale
        return_list.append((new_x,new_y))
        # return_list.append((x3+x1,y3+y1))
    return return_list

def translate_positions(fixed_positions, t_x, t_y):
    new_fixed_positions = {}
    for node in fixed_positions:
        x, y = fixed_positions[node]
        new_fixed_positions[node] = (x + t_x, y + t_y)

    return new_fixed_positions

def get_next_node(node1, circular_nodes):
    for i, node2 in enumerate(circular_nodes):
        if node1 == node2:
            if i+1 < len(circular_nodes):
                return circular_nodes[i+1]
            else:
                return circular_nodes[0]

def get_source_node(node, edge_list):
    for u, v, w in edge_list:
        # if u == node:
        #     return v
        if v == node:
            return u

def draw_graph(component_list, component_edge_list, inter_comp_edge_color, filename, plt_fig):
    G = nx.DiGraph(directed=True)

    total_fixed_positions = {}
    total_fixed_nodes = []

    # print(str(len(component_list)) + " components")
    for ii, component in enumerate(component_list):
        cycled_node_list, uncycled_node_list, edge_list, color = component

        # Format the names to fit inside a node, and change the list accordingly
        adjusted_node_dict = {}
        for node_id in cycled_node_list+uncycled_node_list:
            adjusted_node_dict[node_id] = get_node_name(node_id)
        temp_list = []
        for node_id in cycled_node_list:
            temp_list.append(adjusted_node_dict[node_id])
        cycled_node_list = temp_list
        temp_list = []
        for node_id in uncycled_node_list:
            temp_list.append(adjusted_node_dict[node_id])
        uncycled_node_list = temp_list
        temp_list = []
        for u, v, w in edge_list:
            temp_list.append((adjusted_node_dict[u], adjusted_node_dict[v], w))
        edge_list = temp_list

        x_min, y_min, x_max, y_max = get_border_region(total_fixed_positions)
        # print(round(x_min, 2), round(y_min, 2), round(x_max, 2), round(y_max, 2))
        ww = x_max - x_min
        hh = y_max - y_min
        if ww > hh:
            t_x = 0.0
            t_y = hh
        else:
            t_x = ww
            t_y = 0.0

        circular_positions = get_coordinates_in_circle(len(cycled_node_list))

        G.add_weighted_edges_from(edge_list)

        circular_nodes = cycled_node_list
        fixed_positions = {}
        fixed_nodes = []

        # assign fixed coordinates to circular nodes
        for i, node in enumerate(cycled_node_list):
            fixed_positions[node] = circular_positions[i]
            fixed_nodes.append(node)

        outgoing_nodes_from_circular_points = {}
        for u, v, w in edge_list:
            if u in circular_nodes and v not in circular_nodes:
                if u not in outgoing_nodes_from_circular_points:
                    outgoing_nodes_from_circular_points[u] = []
                outgoing_nodes_from_circular_points[u].append(v)
            elif u not in circular_nodes and v in circular_nodes:
                if v not in outgoing_nodes_from_circular_points:
                    outgoing_nodes_from_circular_points[v] = []
                outgoing_nodes_from_circular_points[v].append(u)

        circular_node_count = len(circular_nodes)
        outside_angle = 2 * math.pi - (circular_node_count - 2) * math.pi / circular_node_count
        # sys.exit()
        for a_node_in_circle in outgoing_nodes_from_circular_points:
            points = get_n_fixed_points_next_to_circle(len(outgoing_nodes_from_circular_points[a_node_in_circle]), outside_angle, fixed_positions[a_node_in_circle], fixed_positions[get_next_node(a_node_in_circle, circular_nodes)])
            for i, node in enumerate(outgoing_nodes_from_circular_points[a_node_in_circle]):
                fixed_positions[node] = points[i]
                fixed_nodes.append(node)

        while len(fixed_nodes) < len(cycled_node_list) + len(uncycled_node_list):
            # print("current_fixed")
            # print(fixed_nodes)
            adj_list = {}
            # for u, v, w in uncycled_edge_list:
            for u, v, w in edge_list:
                if u in fixed_positions and v not in fixed_positions:
                    if u not in adj_list:
                        adj_list[u] = []
                    adj_list[u].append(v)
                elif u not in fixed_positions and v in fixed_positions:
                    if v not in adj_list:
                        adj_list[v] = []
                    adj_list[v].append(u)
            # print "test1"
            for a_node in adj_list:
                # points = get_n_fixed_points_far_from_circle(len(adj_list[a_node]), fixed_positions[a_node], fixed_positions[get_source_node(a_node, cycled_edge_list+uncycled_edge_list)])
                points = get_n_fixed_points_far_from_circle(len(adj_list[a_node]), fixed_positions[a_node], fixed_positions[get_source_node(a_node, edge_list)])
                for i, node in enumerate(adj_list[a_node]):
                    fixed_positions[node] = points[i]
                    fixed_nodes.append(node)
        # print "test2"
        fixed_positions = translate_positions(fixed_positions, t_x, t_y)
        # print "test3"
        for node in fixed_positions:
            total_fixed_positions[node] = fixed_positions[node]
        # print "test4"
        for node in fixed_nodes:
            total_fixed_nodes.append(node)
        # print("component " + str(ii+1))
        # print("fixed_nodes (" + str(len(fixed_nodes)) + ")")
        # print(fixed_nodes)
        # print fixed_positions.keys()

        # edge_labels=dict([((u,v,),d['weight']) for u,v,d in G.edges(data=True)])
        # edge_colors = [color for edge in G.edges()]

        # print G.edges()
        # print "\n" + filename
        # print fixed_nodes
        # print color
        # print edge_colors
        # if ii == 3:
        #     break

    edge_labels = dict([((u,v,),d['weight']) for u,v,d in G.edges(data=True)])
    edge_colors = [color for edge in G.edges()]

    label_options = {
    'font_size': 9, # 10
    'pos':nx.spring_layout(G, pos=total_fixed_positions, fixed=total_fixed_nodes),
    'edge_labels': edge_labels,
    }

    graph_options = {
    'node_color': color,
    'edge_color': edge_colors,
    # 'node_size': [len(v) * 300 for v in G.nodes()],
    'node_size': 1200,
    'node_shape': 'o', # s, o
    'font_size': 9, # 10
    'width': 1,
    # 'arrows': False,
    'arrowstyle': '->',
    'arrowsize': 12,
    'pos':nx.spring_layout(G, pos=total_fixed_positions, fixed=total_fixed_nodes),
    'with_labels': True,    #node label
    # 'edge_labels': edge_labels,
    }
    nx.draw_networkx_edge_labels(G, **label_options)
    # nx.draw_networkx(G, **options)
    nx.draw(G, **graph_options)

    # val_map = {'A': 1.0, 'D': 0.5714285714285714, 'H': 0.0}

    # values = [val_map.get(node, 0.25) for node in G.nodes()]

    # # Specify the edges you want here
    # red_edges = [('A', 'C'), ('E', 'C')]
    # edge_colours = ['black' if not edge in red_edges else 'red' for edge in G.edges()]
    # black_edges = [edge for edge in G.edges() if edge not in red_edges]

    # # Need to create a layout when doing
    # # separate calls to draw nodes and edges
    # pos = nx.spring_layout(G)
    # nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), node_color = values, node_size = 500)
    # nx.draw_networkx_labels(G, pos)
    # nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r')
    # nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=False)
    # plt.show()
    # fig = plt.figure()
    plt.savefig(filename, format="PNG")
    # plt.savefig(filename, dpi=fig.dpi, format="PNG")
    # plt.savefig(filename, dpi=200, format="PNG")
    # plt.clf()
    # plt.close()
    plt_fig.clf()

def get_motif_family_fullname(family_name):
    if family_name.lower() in known_motif_fullname:
        return known_motif_fullname[family_name.lower()]
    else:
        return family_name

def get_loop_source(loop, pdb_organism_details):
    if loop == None or len(loop) == 0 or ":" not in loop:
        return ""
    
    pdb_chain = loop.strip().split(":")[0]
    cluster_source = ""
    if show_cluster_source == True:
        loop = strToNode(loop)
        if loop in loop_cluster_source:
            cluster_source = loop_cluster_source[loop] + ", "
    if pdb_chain in pdb_organism_details:
        org_type = pdb_organism_details[pdb_chain][3][:10] + '.' if len(pdb_organism_details[pdb_chain][3]) > 10 else pdb_organism_details[pdb_chain][3]
        # rna_type = pdb_organism_details[pdb_chain][0][:15] + '...' if len(pdb_organism_details[pdb_chain][0]) > 15 else pdb_organism_details[pdb_chain][0]
        rna_type = pdb_organism_details[pdb_chain][0]

        if len(org_type.strip()) == 0:
            org_type = 'N/A'
        return cluster_source + org_type + ', ' + rna_type
    else:
        return cluster_source + "N/A"

def get_cluster_code(image_name, out):
    # print(image_name)
    pieces = image_name.strip().split("__")
    family_name = pieces[0]
    image_id = pieces[1].strip().split("_")[0]
    if re.match("v\d",family_name) != None:
        family_name = pieces[1]
        image_id = pieces[2].strip().split("_")[0]

    if ('step4' in out or 'step2.0' in out) and image_id != 'all':
        return get_motif_family_short_code(family_name) + "-Sub" + image_id
    return get_motif_family_short_code(family_name) + "_" + image_id

def get_image_and_captions(image_list, show_pdb_info, pdb_organism_details, out):
    new_image_list = []
    title = ""
    for (loop, png, loop_count) in image_list:

        cluster_source = get_loop_source(loop, pdb_organism_details)
        image_name = ".".join(os.path.basename(png).strip().split(".")[:-1])
        title = get_motif_family_fullname(image_name.strip().split("_step")[0])
        if re.match("v\d",title) != None:
            title = get_motif_family_fullname(image_name.strip().split("_")[1])
        cluster_code = get_cluster_code(image_name, out)
        caption = cluster_code if len(cluster_code) < 10 else cluster_code[:9] + '.'

        if output_env == 'local' and loop_count == 1 and allow_two_line_caption == True:
            loop_pdb_ind = convert_a_loop_from_FASTA_to_PDB(loop)
            caption += ';' + loop_pdb_ind + ';' + cluster_source
        else:
            if show_pdb_info == True:
                caption += "," + loop
            else:    
                caption += " ("
                if loop_count > 1:
                    caption += str(loop_count) + " loops"
                else:
                    caption += cluster_source
                caption += ")"
        new_image_list.append((caption, png, loop_count))

    return title, new_image_list

def get_initial_fontsize(width):
    default_string = "ABCDEFGHJKLMNOPQRSTU"
    f_size = default_fontsize
    font_name = "Arial"

    font = ImageFont.truetype(os.path.join(fonts_dir, font_name + ".ttf"), f_size)
    w, h = font.getsize(default_string)
    while w < width and f_size < max_fontsize:
        f_size += 1
        font = ImageFont.truetype(os.path.join(fonts_dir, font_name + ".ttf"), f_size)
        w, h = font.getsize(default_string)

    return f_size

def get_image_caption_size(image_caption, width, initial_fontsize, is_title=False):

    f_size = initial_fontsize
    font_name = "Arial"
    if is_title == True:
        f_size += 5
        font_name = "arialbd"

    font = ImageFont.truetype(os.path.join(fonts_dir, font_name + ".ttf"), f_size)
    w, h = font.getsize(image_caption)
    while w > width and f_size > initial_fontsize - 1:
        f_size -= 1
        font = ImageFont.truetype(os.path.join(fonts_dir, font_name + ".ttf"), f_size)
        w, h = font.getsize(image_caption)
    is_trimmed = False
    while w > width:
        image_caption = image_caption[:-1]
        is_trimmed = True
        w, h = font.getsize(image_caption)
    if is_trimmed == True:
        image_caption = image_caption[:-3] + "..."

    return w, h, image_caption, font

def find_white_boundary_index(pixels, ind_list1, ind_list2, x_first = True):
    boundary_index = 0
    if len(ind_list1) > 0:
        boundary_index = ind_list1[0]
    for i in ind_list1:
        all_white_pixel = True
        for j in ind_list2:         
            if x_first == False:
                r, g, b, a = pixels[j, i]
            else:
                r, g, b, a = pixels[i, j]
            if r != 255 or g != 255 or b != 255:
                all_white_pixel = False
                break
        if all_white_pixel == False:
            break

        boundary_index = i

    return boundary_index

dollar = 1
def crop_image(img, do_save=False):

    if do_save == True:
        output_dir = os.path.join(superimposition_output_dir, 'cropped')

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    pixels = img.load()
    x, y = img.size

    start_x = find_white_boundary_index(pixels, range(x), range(y))
    start_y = find_white_boundary_index(pixels, range(y), range(x), False)
    end_x = find_white_boundary_index(pixels, range(x)[::-1], range(y)[::-1])
    end_y = find_white_boundary_index(pixels, range(y)[::-1], range(x)[::-1], False)

    if end_x == 0:
        end_x = x
    if end_y == 0:
        end_y = y

    # print((start_x, start_y, end_x, end_y))
    img = img.crop((start_x, start_y, end_x, end_y))
    if do_save == True:
        global dollar
        img.save(os.path.join(output_dir, str(dollar) + '.png'))
        dollar += 1
    return img

def is_all_images_generated(image_list):
    for (loop, png, loop_count) in image_list:
        if not os.path.isfile(png):
            return False
    return True

# def is_all_images_generated_name_only(image_list):
#     for png in image_list:
#         if not os.path.isfile(png):
#             return False
#     return True

def create_combined_subfamily_collage(subfamily_dir, image_list, cluster_id, draw_figures, suffix = '', padding_dimension = 40):

    # print image_list
    # sys.exit()

    # if create_subfamily_images == False:
    if draw_figures == False:
        return
    
    # if len(image_list) < 3:
    #     return


    if not os.path.exists(subfamily_dir):
        os.makedirs(subfamily_dir)

    if wait_for_certain_files_to_be_generated(image_list) == False:
        logger.error('All images not generated within max wait time. Collage not created.')
        return

    # wait_for_certain_time_according_to_wait_factor(len(image_list))

    # max_wait_time = 10.0    # seconds
    # wait_time = 0.0
    # while not is_all_images_generated_name_only(image_list):
    #     if wait_time > max_wait_time:
    #         logger.error('All images not generated within max wait time. Collage not created.')
    #         return
    #     logger.info('Waiting for images to be generated.')
    #     time.sleep(.200)
    #     wait_time += .200

    # time.sleep(.100)


    out = os.path.join(subfamily_dir, str(cluster_id) + '_' + suffix + '.png')

    img = Image.open(image_list[-1])
    img = crop_image(img)
    dimension_thumb_x, dimension_thumb_y = img.size
    dimension_thumb_y = dimension_thumb_x

    dimension_x = dimension_thumb_x
    dimension_y = dimension_thumb_y * 3

    img.close()

    target_img = Image.new("RGBA", (int(dimension_x), int(dimension_y)), color=(255,255,255,255))
 
    d = ImageDraw.Draw(target_img)

    running_y = 0

    for row, png in enumerate(image_list):
        img = Image.open(png)
        img.thumbnail((dimension_thumb_x, dimension_thumb_y))
        cur_x, cur_y = img.size
        image_x = int((dimension_thumb_x - cur_x) / 2.0)
        image_y = running_y
        target_img.paste(img, (int(image_x), int(image_y)), mask=img)
        running_y += cur_y + padding_dimension
        img.close()
 
    target_img = crop_image(target_img)
    target_img.save(out)

    wait_for_certain_files_to_be_generated([out], False)


def create_collage(image_list, pdb_organism_details, out, show_pdb_info, is_graph_image, draw_text = True, text_padding = 20, dimension_thumb = 400, padding_dimension = 20, title_height = 0):

    # max_wait_time = 10.0    # seconds
    # wait_time = 0.0
    file_list = []
    for (loop, png, loop_count) in image_list:
        file_list.append(png)

    if wait_for_certain_files_to_be_generated(file_list) == False:
        logger.error('All images not generated within max wait time. Collage not created.')
        return

    # while not is_all_images_generated(image_list):
    #     if wait_time > max_wait_time:
    #         logger.error('All images not generated within max wait time. Collage not created.')
    #         return
    #     logger.info('Waiting for images to be generated.')
    #     time.sleep(.200)
    #     wait_time += .200

    # time.sleep(.100)

    title, image_list = get_image_and_captions(image_list, show_pdb_info, pdb_organism_details, out)
    if is_graph_image == True:
        dimension_thumb = 800
    # You can add a check to guarantee that the list of images is not bigger than 64...
    no_of_files = len(image_list)
    root = math.sqrt(no_of_files)
    # print root

    if collage_shape == 'rectangle':
        root /= 2.0

    root = math.ceil(root)

    if len(image_list) > 90:
        dimension_thumb = (dimension_thumb) / 4   #100    150   200
    elif len(image_list) > 56:
        dimension_thumb = (dimension_thumb * 3) / 8   #150  225 300
    elif len(image_list) > 30:
        dimension_thumb = (dimension_thumb) / 2   #200    300    400
    elif len(image_list) > 12:
        dimension_thumb = (dimension_thumb * 3) / 4   #300  450    600

    initial_fontsize = get_initial_fontsize(dimension_thumb)
    # factor = int((root - 1)/4) + 1
    # dimension thumb size changes for more than 4+, 8+, 12+, etc rows
    # dimension_thumb = int(dimension_thumb / factor)

    count_x = int(math.ceil(no_of_files / float(root)))
    dimension_x = count_x * dimension_thumb + padding_dimension * (count_x + 1)
    dimension_y = int(root) * dimension_thumb  + padding_dimension * (int(root) + 1) + title_height
    if draw_text == True:
        if output_env == 'local' and len(image_list) == sum(map(lambda x: x[2], image_list)) and allow_two_line_caption == True:
            dimension_y += 2 * text_padding * int(root)
        else:
            dimension_y += text_padding * int(root)

    target_img = Image.new("RGBA", (int(dimension_x), int(dimension_y)), color=(255,255,255,255))
 
    d = ImageDraw.Draw(target_img)

    if title_height > 0:
        w, h, title, font = get_image_caption_size(title, dimension_x, initial_fontsize, True)
        text_x = int((dimension_x - w) / 2.0)
        text_y = int((title_height - h) / 2.0)
        d.text((text_x, text_y), title, fill=(0, 0, 0), font=font)

    for k, (image_caption, png, loop_count) in enumerate(image_list):
        image_basename = os.path.basename(png)
        row, col = divmod(k, int(count_x))
        img = Image.open(png)

        img = crop_image(img)

        img.thumbnail((dimension_thumb, dimension_thumb))
        cur_x, cur_y = img.size
        image_x = dimension_thumb*col + padding_dimension * (col + 1) + int((dimension_thumb - cur_x) / 2.0)
        image_y = dimension_thumb*row + padding_dimension * (row + 1) + int((dimension_thumb - cur_y) / 2.0) + title_height
        if draw_text == True:
            if image_caption.count(';') == 2 and loop_count == 1 and allow_two_line_caption == True:
                image_y += 2 * text_padding * row
            else:
                image_y += text_padding * row

        target_img.paste(img, (int(image_x), int(image_y)), mask=img)
        img.close()
        
        if draw_text == True:
            text_x = int(dimension_thumb*col + padding_dimension * (col + 1.1))
            text_y = int(dimension_thumb*(row + 1) + padding_dimension * (row + 1.1) + text_padding * row) + title_height
            d = ImageDraw.Draw(target_img)

            if image_caption.count(';') == 2 and loop_count == 1 and allow_two_line_caption == True:
                text_y += int(text_padding * row)

                caption_pieces = image_caption.strip().split(';')
                caption_line1 = ', '.join(caption_pieces[:2])
                caption_line2 = caption_pieces[2]

                w, h, caption_line1, font = get_image_caption_size(caption_line1, dimension_thumb, initial_fontsize)
                text_x1 = text_x + int((dimension_thumb - w) / 2.0)
                text_y1 = text_y + int((text_padding - h) / 2.0)
                d.text((text_x1, text_y1), caption_line1, fill=(0, 0, 0), font=font)

                text_y1 += int(h * 1.1)
                w, h, caption_line2, font = get_image_caption_size(caption_line2, dimension_thumb, initial_fontsize)
                text_x2 = text_x + int((dimension_thumb - w) / 2.0)
                text_y2 = text_y1 + int((text_padding - h) / 2.0)# + int(text_padding * 1.5)
                d.text((text_x2, text_y2), caption_line2, fill=(0, 0, 0), font=font)
            else:
                w, h, image_caption, font = get_image_caption_size(image_caption, dimension_thumb, initial_fontsize)
                text_x += int((dimension_thumb - w) / 2.0)
                text_y += int((text_padding - h) / 2.0)
                d.text((text_x, text_y), image_caption, fill=(0, 0, 0), font=font)
 
    target_img.save(out)

    wait_for_certain_files_to_be_generated([out], False)

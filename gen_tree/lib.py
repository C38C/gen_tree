import scipy.spatial
import scipy.special
import math
import time
import random
import numpy as np
import os


def dist_2d(point1, point2):
    return ((point1[0] - point2[0]) ** 2.0 + (point1[1] - point2[1]) ** 2.0 + \
        (point1[2]-point2[2])**2.0)**0.5


def area_triangle(point1, point2, point3):  # using heron's formula
    a = dist_2d(point1, point2)
    b = dist_2d(point2, point3)
    c = dist_2d(point3, point1)
    s = (a + b + c) / 2.0
    return (s * (s-a) * (s-b) * (s-c)) ** 0.5


def centroid_triangle(point1, point2, point3):
    return ((point1[0] + point2[0] + point3[0]) / 3.0, \
        (point1[1] + point2[1] + point3[1]) / 3.0, \
        (point1[2] + point2[2] + point3[2]) / 3.0)


def area_hemisphere(width):
    return 2 * math.pi * (width/2.0) ** 2.0


def area_ellipsoid(c, a):
    e = (1 - (c**2/a**2))**0.5
    return 2 * math.pi * a**2 * (1 + (1 - e**2) / e * math.atanh(e))


def area_prolate_ellipsoid(c, a):
    e = (1 - (a ** 2 / c ** 2)) ** 0.5
    return (2 * math.pi * a ** 2 * (1 + (c / (a * e) * math.asin(e))))


def area_cone(r, h):
    return r * math.pi * ((h ** 2 + r ** 2) ** 0.5)

# Generate vertices using parametric equation for cone


def generate_cone_vertices(vertices_count, height, radius):
    vertices = []
    while len(vertices) < vertices_count:
        # create z points with a lesser probability of occurence with increasing height.
        u = np.random.triangular(0, 0, height)
        theta = random.random() * math.pi * 2

        x = ((height - u) / height) * radius * math.cos(theta)
        y = ((height - u) / height) * radius * math.sin(theta)
        z = u

        vertices.append((x, y, z))
    return vertices


def generate_ellipsoid_vertices(vertices_count, height, radius):
    vertices = []
    while len(vertices) < vertices_count:
        x1 = random.uniform(-1.0, 1.0)
        x2 = random.uniform(-1.0, 1.0)

        t1 = (x1 ** 2.0 + x2 ** 2.0)
        t2 = (1.0 - 2.0 * (t1))

        if t1 < 1.0 and t2 >= 0.0:
            x = 2.0 * x1 * (1.0 - x1 ** 2.0 - x2 ** 2.0) ** 0.5
            y = 2.0 * x2 * (1.0 - x1 ** 2.0 - x2 ** 2.0) ** 0.5
            z = t2
            mu = np.linalg.norm(
                [radius*height*y, radius*radius*z, radius*height*x])
            mumax = max([radius*height, radius*radius])

            if random.random() < mu/mumax:
                vertices.append((radius*x, radius*y, height*z))
    return vertices


def generate_hemisphere_vertices(vertices_count, width):
    vertices = []  # list containing all vertices
    scale_factor = width/2.0

    while len(vertices) < vertices_count:  # generate vertices using Marsaglia's method
        x1 = random.uniform(-1.0, 1.0)  # random floating point within range
        x2 = random.uniform(-1.0, 1.0)

        t1 = (x1**2.0 + x2**2.0)
        t2 = (1.0-2.0*(t1))

        if t1 < 1.0 and t2 >= 0.0:
            x = 2.0*x1*(1.0-x1**2.0-x2**2.0)**0.5
            y = 2.0*x2*(1.0-x1**2.0-x2**2.0)**0.5
            z = t2

            vertices.append((scale_factor*x, scale_factor*y, scale_factor*z))
    return vertices


def triangulate(vertices):
    v2d = [(v[0], v[1]) for v in vertices]
    return scipy.spatial.Delaunay(v2d)


def triangulate_convex_hull(vertices):
    return scipy.spatial.ConvexHull(vertices, qhull_options="Qt")


def write_point_cloud(filename, vertices):
    if not os.path.exists("output"):
        os.makedirs("output")
    with open("output/" + filename + ".csv", 'w') as f:
        f.write('x,y,z\n')
        for point in vertices:
            f.write('%.8f,%.8f,%.8f\n' % point)


def write_obj_file(filename, vertices, faces):
    if not os.path.exists("output"):
        os.makedirs("output")
    with open("output/" + filename + ".obj", 'w') as f:
        f.write('#gen_tree OBJ output...\n')
        for v in vertices:
            # for v in v2d:
            f.write('v %s %s %s\n' % (v[0], v[1], v[2]))
            # f.write('v %s %s\n' % tuple(v))
        f.write('g mesh\n')

        for face in faces:
            f.write('f %i// %i// %i//\n' % (face[0]+1, face[1]+1, face[2]+1))
# Write


def write_rad_file(filename, vertices, faces):
    if not os.path.exists("output"):
        os.makedirs("output")
    with open("output/" + filename + ".rad", 'w') as f:
        for (n, face) in enumerate(faces):
            f.write("tree_material polygon %s.%i \n" % (filename, n))
            f.write("0 \n0 \n9 ")
            for vtx_index in face:
                f.write("\t\t%.8f \t\t%.8f \t\t%.8f\n" % vertices[vtx_index])
            f.write("\n")


def plotting_point_cloud(vertices):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    xs = [v[0] for v in vertices]
    ys = [v[1] for v in vertices]
    zs = [v[2] for v in vertices]
    ax.scatter(xs, ys, zs)

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()

# get_good_triangles removes the bottom and base tringles of the geometry


def get_good_triangles(vertices, faces, h, max_length):
    triangle_faces = []
    for face in faces:
        good_triangle = True
        # Checking if any of the triangle faces lie close to the base of the cone.
        # To be precise at a height of 1/50th (approx) of the radius of the cone.
        if vertices[face[0]][2] + vertices[face[1]][2] + vertices[face[2]][2] < (0.025 * h):
            good_triangle = False

        else:
            # larger vertex count, smaller the triangles and smaller vertex count larger the triangles.
            for face_index in range(len(face)):
                # Checking the distance between any one side of the triangular face.
                # Change this max_length factor if you smaller or bigger triangles.
                if dist_2d(vertices[face[face_index]], vertices[face[face_index-1]]) > max_length:
                    good_triangle = False
                    break

        if good_triangle:
            triangle_faces.append(face.tolist())
    return triangle_faces


def get_triangle_areas(vertices, faces):
    return [area_triangle(vertices[face[0]], vertices[face[1]], vertices[face[2]]) for face in faces]


def get_triangle_centroids(vertices, faces):
    return [centroid_triangle(vertices[face[0]], vertices[face[1]], vertices[face[2]]) for face in faces]


def select_faces(vertices, faces, area, gap_percentage, cluster_density=3.0):
    cluster_distance_factor = (area / len(faces)) ** 0.5
    face_keep_indices = []  # indices of faces to keep
    triangle_centroids = get_triangle_centroids(vertices, faces)
    triangle_areas = get_triangle_areas(vertices, faces)

    # compute KDTree of point distances
    kd_tree = scipy.spatial.KDTree(triangle_centroids)

    goal_to_keep = (1 - (gap_percentage ** 0.5)) * area  # area of tree to be kept
    total_kept = 0.0  # tracks kept area

    print("goal_to_keep %f" % goal_to_keep)

    start_time = time.time()
    total_time = 0
    print(start_time)
    print("starting the loop")

    while total_kept < goal_to_keep and (total_time < 50000):
        total_time += 1
        face_index = random.randint(0, len(faces)-1)  # random face on tree
        # random selection of fill option: single leaf vs. cluster
        choice = random.randint(0, 1)

        if choice == 0:  # single, lonely leaf
            if face_index not in face_keep_indices:  # if face isn't already filled
                face_keep_indices.append(face_index)
                total_kept += triangle_areas[face_index]
        elif choice == 1:  # cluster of leaves
            # cluster distance
            distance = random.random() * cluster_density * cluster_distance_factor
            # increasing this factor from 5.0 to 10.0 makes the loosely filled cluster -- that is trees that are not dense
            # descreasing the number to 2.0 gives an almost uniformly distributed look without much clusters
            # close faces with centroids within cluster distance
            close_faces = kd_tree.query_ball_point(
                triangle_centroids[face_index], distance)
            for cf in close_faces:  # for each close face
                # 2/7 chance of hole in cluster of leaves
                if (cf not in face_keep_indices) and (random.randint(0, 8) > 1):
                    face_keep_indices.append(cf)
                    total_kept += triangle_areas[cf]

    print("total_kept: %f" % total_kept)
    print("Time elapsed")
    return face_keep_indices


def main_run(vertices, gap_percentage, max_length, shape_name, height, ideal_area,
             cluster_density=3.0, convex_hull=False):
    # Write the point cloud generated into a csv file.
    write_point_cloud(shape_name + "PointCloud", vertices)

    # Create triangulated mesh using Delaunay
    tmp_triangles = []
    if convex_hull:
        tmp_triangles = triangulate_convex_hull(vertices)
    else:
        tmp_triangles = triangulate(vertices)

    # Getting faces and vertices of the mesh
    tmp_triangle_faces = tmp_triangles.simplices
    tmp_triangle_vertices = tmp_triangles.vertices

    # Keeping the good triangles from the mesh.
    triangle_faces = get_good_triangles(
        vertices, tmp_triangle_faces, height, max_length)
    print("after removing count = %d" % len(triangle_faces))

    # Getting the final triangulated geometry in accordance with gap percentage.
    # Increase cluster_density for denser cluster. Decrease cluster_density for even looking crown.
    face_keep_indices = select_faces(
        vertices, triangle_faces, ideal_area, gap_percentage,
        cluster_density=cluster_density)

    # =========================================
    # output final mesh as .obj and .rad files
    # =========================================
    kept_faces = [triangle_faces[i] for i in face_keep_indices]
    write_obj_file(shape_name, vertices, kept_faces)
    write_rad_file(shape_name, vertices, kept_faces)


def gen_hemisphere(file_name="Hemisphere", vertices_count=20000, gap_percentage=0.136, width=5.0):
    width = float(width)
    vertices = generate_hemisphere_vertices(vertices_count, width)

    # max_length is based on the 1/2 circumferance length divided by the sqrt of the vertex count
    max_length = 4.0 * math.pi * (width/2.0) / (vertices_count**0.5)

    ideal_area_hemisphere = area_hemisphere(width)
    print ("area %f" % ideal_area_hemisphere)

    main_run(vertices, gap_percentage, max_length, file_name, width,
             ideal_area_hemisphere, cluster_density=2.0, convex_hull=True)


def gen_ellipsoid_prolate(file_name="EllProlate", vertices_count=20000, gap_percentage=0.059, height=7.8,
                          radius=5.35):
    height = float(height)
    radius = float(radius)
    vertices = generate_ellipsoid_vertices(vertices_count, height, radius)

    # Area of half an prolate ellipsoid.
    ideal_area_ellipsoid = area_prolate_ellipsoid(height, radius) / 2

    # now remove the bottom triangles of the triangulated hemisphere
    max_length = 4 * math.pi * radius / (vertices_count**0.5)
    # max_length is based on the 1/2 circumferance length divided by the sqrt of the vertex count

    main_run(vertices, gap_percentage, max_length, file_name, height,
             ideal_area_ellipsoid, cluster_density=7.0, convex_hull=True)


def gen_ellipsoid_oblate(file_name="EllOblate", vertices_count=20000, gap_percentage=0.03, height=5.2,
                         radius=5.7):
    height = float(height)
    radius = float(radius)
    vertices = generate_ellipsoid_vertices(vertices_count, height, radius)

    # Area of half an oblate ellipsoid.
    ideal_area_ellipsoid = area_ellipsoid(height, radius) / 2

    # Remove the bottom triangles of the triangulated hemisphere
    max_length = 4.0 * ((math.pi * radius) / (vertices_count**0.5))
    # Notes:
    # max_length is based on the 1/2 circumferance length divided by the sqrt of the vertex count
    # the smaller this factor the even looking triangles make the mesh
    # the smaller this number the difficult it is to be a good triangle.

    main_run(vertices, gap_percentage, max_length, file_name, height,
             ideal_area_ellipsoid, cluster_density=7.0, convex_hull=False)


def get_cone(file_name="Cone", vertices_count=10000, gap_percentage=0.18, height=8.5, radius=2.5):
    height = float(height)
    radius = float(radius)
    vertices = generate_cone_vertices(vertices_count, height, radius)

    # Remove the bottom triangles of the triangulated hemisphere
    max_length = 6.0 * math.pi * radius / (vertices_count ** 0.5)

    ideal_area_cone = area_cone(radius, height)
    print("area: %f" %ideal_area_cone)

    main_run(vertices, gap_percentage, max_length, file_name,
             height, ideal_area_cone, cluster_density=3.0)

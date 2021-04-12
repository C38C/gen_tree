import lib

#To generate a hemispherical tree crown. Use this shape when the crown's height and radius are similar to each other.
lib.gen_hemisphere(file_name="Hemisphere", vertices_count=20000, gap_percentage=0.136, width=5.0)

#To generate a truncated prolate ellipsoid tree crown. Use this shape when the crown's height > radius
lib.gen_ellipsoid_prolate(file_name="EllProlate", vertices_count=20000, gap_percentage=0.059, height=7.8, radius=5.35)

#To generate a truncated oblate ellipsoid tree crown. Use this shape when the crown's height < radius
lib.gen_ellipsoid_oblate(file_name="EllOblate", vertices_count=20000, gap_percentage=0.07, height=8.0, radius=10.0)

#To generate a conical tree crown. Use this shape when the tree crown is broad at the base and narrows to the top.
lib.get_cone(file_name="Cone", vertices_count=10000, gap_percentage=0.18, height=8.5, radius=2.5)

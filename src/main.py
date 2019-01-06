import lib

# #To generate a hemispherical tree crown. Use this shape when the crown's height and radius are similar to each other.
lib.genHemisphere(vertices_count=20000, gapPercentage=0.136, width=5.0)
#
# #To generate a truncated prolate ellipsoid tree crown. Use this shape when the crown's height > radius
lib.genEllipsoidProlate(vertices_count=20000, gapPercentage=0.059, height=7.8, radius=5.35)

#To generate a truncated oblate ellipsoid tree crown. Use this shape when the crown's height < radius
lib.genEllipsoidOblate(vertices_count=20000, gapPercentage=0.07, height=8.0, radius=10.0)

#To generate a conical tree crown. Use this shape when the tree crown is broad at the base and narrows to the top.
lib.genCone(vertices_count=10000, gapPercentage=0.18, height=8.5, radius=2.5)

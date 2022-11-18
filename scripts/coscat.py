import shapefile
from dbfread import DBF

# shapes = shapefile.Reader("Continental_Shelf.shp")
db = DBF("Continental_Shelf.dbf")

for record in DBF("Continental_Shelf.dbf"):   
    a = record
    pass

pass
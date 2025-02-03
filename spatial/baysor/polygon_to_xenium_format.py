import json

def reformat_seg_polygons(data):
"""
Reformat 2d-polygon file to be imported into Xenium Explorer 
@param data: segmentation_polygons_2d.json Baysor output file read into a dictionary
@returns processed: data dictionary reformatted according to schema below
"""
    geometries = []
    for row in data.get('features'):
        geometries_dict = row.get('geometry')
        cell_id = row.get('id')
        cell_id = cell_id.split("-")[1]
        geometries_dict.update({'cell': cell_id})
        geometries.append(geometries_dict)

    processed = {'geometries': geometries}
    processed.update({'type': 'GeometryCollection'})
    return processed                  


"""
This is the format we need the data in for importing into xenium:
{
  "geometries": [
    {
      "coordinates": [
        [
          [15.891722, 21.441122],
          [13.95383, 21.565407],
          [14.233208, 25.427168],
          [15.927658, 26.113512],
          [17.464699, 25.211208],
          [16.642977, 22.83285],
          [15.891722, 21.441122]
        ]
      ],
      "type": "Polygon",
      "cell": 5
    }
  ],
  "type": "GeometryCollection"
}
"""


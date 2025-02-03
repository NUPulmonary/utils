"""
Update Baysor 2d-polygon file so that it can be imported into Xenium Explorer
"""

import json

with open('/projects/b1038/Pulmonary/lcusick/projects/baysor/mlung-test/output/XETG00190_2/segmentation_polygons_2d.json', 'r') as f:
    data = json.load(f)

geometries = []
for row in data.get('features'):
    geometries_dict = row.get('geometry')
    cell_id = row.get('id')
    cell_id = cell_id.split("-")[1]
    geometries_dict.update({'cell': cell_id})
    geometries.append(geometries_dict)

processed = {'geometries': geometries}
processed.update({'type': 'GeometryCollection'})
                      
with open('/projects/b1038/Pulmonary/lcusick/projects/baysor/mlung-test/output/XETG00190_2/segmentation_polygons_2d_processed.json', 'w') as f:
    json.dump(processed, f)

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


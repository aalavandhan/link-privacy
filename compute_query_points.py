import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans

datasets_parameters = {
  "data_Gowalla" : {
    'eps': 40,
    'min_samples': 5,
  },
  "data_GowallaFull" : {
    'eps':  40,
    'min_samples': 5,
  },
  "data_Shanghai" : {
    'eps': 40,
    'min_samples': 100,
  },
}

#for d in sorted(datasets_parameters.keys()):
for d in ["data_Shanghai"]:
  params = datasets_parameters[d]
  print "Processing {0}".format(d)
  locations = pd.read_csv(d + "/checkins.txt", sep=" ", header=None, index_col=None, names=[ "user_id", "lat", "lon", "location_id", "date", "time"])
  print "Grouping locations"
  coordinates = locations[['location_id', 'lat', 'lon']].groupby(['location_id']).mean()[[ "lat", "lon" ]].reset_index()
  print "{0}\tNumber of locations : {1}".format(d, len(coordinates))

  # clusters = DBSCAN(eps=params['eps'], min_samples=params['min_samples'], leaf_size=10, algorithm="kd_tree").fit(coordinates.as_matrix())
  clusters = KMeans(n_clusters=1200).fit(coordinates.as_matrix())

  lat = [ ]
  lon = [ ]

  labels = pd.Series(clusters.labels_)

  for l in labels.unique():
      cluster_locations = labels[ labels == l ]
      e = np.random.choice( cluster_locations.index )

      r = coordinates.ix[e]
      lat.append(    r['lat']    )
      lon.append(    r['lon']    )

  reduced_locations = pd.DataFrame({ 'lat': lat, 'lon': lon })
  reduced_locations['day'] = 0
  print "{0}\tNumber of clusters : {1}".format(d, len(reduced_locations))
  reduced_locations.to_csv(d + "/queries-initial.txt", header=None, index=None, sep="\t")

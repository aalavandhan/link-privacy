import pandas as pd
from sklearn.cluster import DBSCAN

datasets_parameters = {
  "data_Gowalla" : {
    'eps': 40,
    'min_samples': 5,
  },
  "data_GowallaFull" : {
    'eps': 40,
    'min_samples': 5,
  },
  "data_Shanghai" : {
    'eps': 0.5,
    'min_samples': 5,
  },
}

for d in datasets_parameters.keys():
  params = datasets_parameters[d]

  locations = pd.read_csv(d + "/checkins.txt", sep=" ", header=None, index_col=None, names=[ "user_id", "lat", "lon", "location_id", "date", "time"])
  clusters = DBSCAN(eps=params['eps'], min_samples=params['min_samples'], metric='euclidean', algorithm='auto', leaf_size=30, p=None, n_jobs=1).fit(locations[[ "lat", "lon" ]])

  lat = [ ]
  lon = [ ]
  labels = pd.Series(clusters.labels_)

  for l in labels.unique():
      cluster_locations = labels[ labels == l ]
      e = np.random.choice( cluster_locations.index )
      r = locations.ix[e]
      lat.append(    r['lat']    )
      lon.append(    r['lon']    )

  print "{0}\tNumber of clusters : {1}".format(d, len(reduced_locations))
  reduced_locations = pd.DataFrame({ 'lat': lat, 'lon': lon })
  reduced_locations.to_csv(d + "/queries-initial.txt", header=None, names=[ "lat", "lon" ], index=None, sep="\t")

data_set=0

# Running locality generation
./SPFP $data_set compute-locality 0 0 0 0 0 0 0 > results-$data_set/locality-generation
# Running histograms
./SPFP $data_set compute-histograms 0 2.9 20 0 0 0 0 > results-$data_set/histogram-generation

# No noise
./SPFP $data_set ebm 0 0 0 0 1200 0 0 > results-$data_set/plain-ebm

# Running experiments ( Ideal grouping )
# Ideal spatial grouping
./SPFP $data_set grouped-ebm 0 0 0 3.3 1200 0 0 > results-$data_set/ideal-grouping-sp-3.3
./SPFP $data_set grouped-ebm 0 0 0 10 1200 0 0 > results-$data_set/ideal-grouping-sp-10
./SPFP $data_set grouped-ebm 0 0 0 50 1200 0 0 > results-$data_set/ideal-grouping-sp-50
./SPFP $data_set grouped-ebm 0 0 0 100 1200 0 0 > results-$data_set/ideal-grouping-sp-100
./SPFP $data_set grouped-ebm 0 0 0 200 1200 0 0 > results-$data_set/ideal-grouping-sp-200
./SPFP $data_set grouped-ebm 0 0 0 400 1200 0 0 > results-$data_set/ideal-grouping-sp-400
./SPFP $data_set grouped-ebm 0 0 0 800 1200 0 0 > results-$data_set/ideal-grouping-sp-800

# Ideal time grouping
./SPFP $data_set grouped-ebm 0 0 0 0 600 0 0  > results-$data_set/ideal-grouping-tm-600
./SPFP $data_set grouped-ebm 0 0 0 0 1200 0 0  > results-$data_set/ideal-grouping-tm-1200
./SPFP $data_set grouped-ebm 0 0 0 0 2400 0 0  > results-$data_set/ideal-grouping-tm-2400
./SPFP $data_set grouped-ebm 0 0 0 0 3600 0 0  > results-$data_set/ideal-grouping-tm-3600
./SPFP $data_set grouped-ebm 0 0 0 0 7200 0 0  > results-$data_set/ideal-grouping-tm-7200
./SPFP $data_set grouped-ebm 0 0 0 0 14400 0 0  > results-$data_set/ideal-grouping-tm-14400
./SPFP $data_set grouped-ebm 0 0 0 0 28800 0 0  > results-$data_set/ideal-grouping-tm-28800


# Gaussian spatial Noise
for spatial_noise in 100 200 300 400 500 600 700; do
  for spatial_grouping in 0.33; do
    ./SPFP $data_set gaussian-v-ebm 1 $spatial_noise 0 $(( $spatial_noise * $spatial_grouping )) 1200 1 0 > results-$data_set/gaussian-spatial-$spatial_noise-$spatial_grouping
  done
done


# Gaussian temporal Noise
for temporal_noise in 0 20 40 60 80 100 120 140 160 180; do
  for time_grouping in 0 0.25 0.5 1 1.5; do
    ./SPFP $data_set gaussian-v-ebm 1 0 $(( $temporal_noise * 60 )) 0 $(( $temporal_noise * $time_grouping * 60 )) 2 0 > results-$data_set/gaussian-temporal-$temporal_noise-$time_grouping
  done
done

# Both Spatial and temporal gaussian noise functions ( 1/3 grouping )
./SPFP $data_set gaussian-v-ebm 1 400 36000 133 12000 0 0 > results-$data_set/gaussian-both-400-36000-133-12000

# Only spatial noise function
for spatial_noise in 100 200 300 400 500 600 700; do
  for spatial_grouping in 0 0.33 0.5 1 1.5; do
    ./SPFP $data_set comb-v-ebm 1 $spatial_noise 0 $(( $spatial_noise * $spatial_grouping )) 1200 1 0 > results-$data_set/comb-spatial-$spatial_noise-$spatial_grouping
  done
done

# Only temporal noise function ( Minutes )
for temporal_noise in 0 20 40 60 80 100 120 140 160 180; do
  for time_grouping in 0 0.25 0.5 1 1.5; do
    ./SPFP $data_set comb-v-ebm 1 0 $(( $temporal_noise * 60 )) 0 $(( $temporal_noise * $time_grouping * 60 )) 2 0 > results-$data_set/comb-temporal-$temporal_noise-$time_grouping
  done
done

# Both Spatial and temporal noise functions ( 1/3 grouping )
./SPFP $data_set comb-v-ebm 1 400 3600 133 3600 0 0 > results-$data_set/comb-both-400-3600-133-1200

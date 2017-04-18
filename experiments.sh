# Running locality generation
# ./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt compute-locality 0 0 0 0 0 0 0 > gowalla-full-experiments/locality-generation

# Running histograms
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt compute-histograms 0 2.9 20 0 0 0 0 > gowalla-full-experiments/histogram-generation

# No noise
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt ebm 0 0 0 0 1200 0 0 > gowalla-full-experiments/plain-ebm

# Running experiments ( Ideal grouping )
# Ideal spatial grouping
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 3.3 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-3.3
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 10 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-10
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 50 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-50
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 100 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-100
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 200 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-200
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 400 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-400
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 800 1200 0 0 > gowalla-full-experiments/ideal-grouping-sp-800

# Ideal time grouping
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 600 0 0  > gowalla-full-experiments/ideal-grouping-tm-600
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 1200 0 0  > gowalla-full-experiments/ideal-grouping-tm-1200
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 2400 0 0  > gowalla-full-experiments/ideal-grouping-tm-2400
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 3600 0 0  > gowalla-full-experiments/ideal-grouping-tm-3600
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 7200 0 0  > gowalla-full-experiments/ideal-grouping-tm-7200
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 14400 0 0  > gowalla-full-experiments/ideal-grouping-tm-14400
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt grouped-ebm 0 0 0 0 28800 0 0  > gowalla-full-experiments/ideal-grouping-tm-28800

# Only spatial noise function
for spatial_noise in 100 200 300 400 500 600 700; do
  for spatial_grouping in 0 0.5 1 1.5 2 2.5 3; do
    ./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt comb-v-ebm 1 $spatial_noise 0 $(( $spatial_noise * $spatial_grouping )) 1200 1 0 > gowalla-full-experiments/comb-spatial-$spatial_noise-$spatial_grouping
  done
done

# Only temporal noise function ( Minutes )
for temporal_noise in 0 20 40 60 80 100 120 140 160; do
  for time_grouping in 0 0.5 1 1.5 2 2.5 3; do
    ./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt comb-v-ebm 1 0 $(( $temporal_noise * 60 )) 0 $(( $temporal_noise * $time_grouping * 60 )) 2 0 > gowalla-full-experiments/comb-temporal-$temporal_noise-$time_grouping
  done
done

# Both Spatial and temporal noise functions ( 1/3 grouping )
./SPFP data_GowallaFull/socialGraph.txt data_GowallaFull/checkins.txt queries.txt comb-v-ebm 1 400 360000 133 120000 0 0 > gowalla-full-experiments/comb-both-400-360000-133-120000

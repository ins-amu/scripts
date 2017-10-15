# SCript to create custom connectivity matrix

export PRD=/usr/local/freesurger/subjects/R005484092_2

# export tracts to tmp txt file
mkdir -f $PRD/connectivity/tmp_ascii_tck

tckconvert $PRD/connectivity/whole_brain_post.tck $PRD/connectivity/tmp_ascii_tck/output-[].txt -nthreads 16 


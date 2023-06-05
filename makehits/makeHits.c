#include <stdio.h>

#include <makeHits.h>

int makeHits(itape_header_t *event)
{
  int total=0;
  total+=make_lgd_hits(event, NULL, 1000);
  total+=make_bsd_hits(event, NULL, 1000);
  total+=make_tagger_hits(event, NULL, 1000);
  total+=make_cpv_hits(event, NULL, 1000);
  total+=make_upv_hits(event, NULL, 1000);
  total+=make_bgv_hits(event, NULL, 1000);
  total+=make_time_list(event, NULL, MAX_TIME_LIST_LENGTH);

  suppress_dead_channels(event);
  
  total+=make_bsd_pixels(event, NULL, 1000);
  total+=make_lgd_clusters(event, NULL, 1000);
  total+=make_bgv_clusters(event, NULL, 1000);
  total+=make_recoil(event,NULL, 1000);
  total+=make_photons(event, NULL, 1000);
  total+=make_tagger_clusters(event,NULL,MAX_TAGGER_CLUSTERS);
  total+=make_tagger_photons(event,NULL,MAX_TAGGER_PHOTONS);
  return(total);
}

int suppress_dead_channels(itape_header_t *event)
{
  bsd_hits_t *bsd_hits;
  bgv_hits_t *bgv_hits;
  cpv_hits_t *cpv_hits;
  upv_hits_t *upv_hits;
  lgd_hits_t *lgd_hits;
  tagger_hits_t *tagger_hits;

  if (bsd_hits=data_getGroup(event, GROUP_BSD_HITS, 0)) {
    suppress_dead_bsd_channels(bsd_hits);
  }
  if (bgv_hits=data_getGroup(event, GROUP_BGV_HITS, 0)) {
    suppress_dead_bgv_channels(bgv_hits);
  }
  if (cpv_hits=data_getGroup(event, GROUP_CPV_HITS, 0)) {
    suppress_dead_cpv_channels(cpv_hits);
  }
  if (upv_hits=data_getGroup(event, GROUP_UPV_HITS, 0)) {
    suppress_dead_upv_channels(upv_hits);
  }
  if (lgd_hits=data_getGroup(event, GROUP_LGD_HITS, 0)) {
    suppress_dead_lgd_channels(lgd_hits);
  }
  if (tagger_hits=data_getGroup(event, GROUP_TAGGER_HITS, 0)) {
    suppress_dead_tagger_channels(tagger_hits);
  }
  return 0;
}

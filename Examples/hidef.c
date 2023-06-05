#include <string.h>

#include <hituple.h>
#include <cernlib.h>

void appendnt_bsd_hits(int id, nt_bsd_hits_t *p_bsd)
{
  if (hexist_(&id)) {
    if (p_bsd != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"bsd_hits",p_bsd,NT_FORM_BSD_HIT);
      }
      else {
        hbname(id,"bsd_hits",p_bsd,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_bsd_hits:");
      fprintf(stderr," p_bsd pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_bsd_hits - ntuple id %d does not exist!\n",
            id);
  } 
}

void appendnt_bsd_pixs(int id, nt_bsd_pixs_t *p_pix)
{
  if (hexist_(&id)) {
    if (p_pix != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"bsd_pixs",p_pix,NT_FORM_BSD_PIX);
      }
      else {
        hbname(id,"bsd_pixs",p_pix,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_bsd_pixs:");
      fprintf(stderr," p_pix pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_bsd_pixs - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_bgv_hits(int id, nt_bgv_hits_t *p_bgv)
{
  if (hexist_(&id)) {
    if (p_bgv != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"bgv_hits",p_bgv,NT_FORM_BGV_HIT);
      }
      else {
        hbname(id,"bgv_hits",p_bgv,"$SET");
      }
    } 
    else {
      fprintf(stderr,"WARNING: appendnt_bgv_hits:");
      fprintf(stderr," p_bgv pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_bgv_hits - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_cpv_hits(int id, nt_cpv_hits_t *p_cpv)
{
  if (hexist_(&id)) {
    if (p_cpv != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"cpv_hits",p_cpv,NT_FORM_CPV_HIT);
      }
      else {
        hbname(id,"cpv_hits",p_cpv,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_cpv_hits:");
      fprintf(stderr," p_cpv pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_cpv_hits - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_upv_hits(int id, nt_upv_hits_t *p_upv)
{
  if (hexist_(&id)) {
    if (p_upv != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"upv_hits",p_upv,NT_FORM_UPV_HIT);
      }
      else {
        hbname(id,"upv_hits",p_upv,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_upv_hits:");
      fprintf(stderr," p_upv pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_upv_hits - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_tag_hits(int id, nt_tag_hits_t *p_tag)
{

  if (hexist_(&id)) {
    if (p_tag != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"tag_hits",p_tag,NT_FORM_TAG_HIT);
      }
      else {
        hbname(id,"tag_hits",p_tag,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_tag_hits:");
      fprintf(stderr," p_tag pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_tag_hits - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_lgd_hits(int id, nt_lgd_hits_t *p_lgd)
{
  if (hexist_(&id)) {
    if (p_lgd != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"lgd_hits",p_lgd,NT_FORM_LGD_HIT);
      }
      else {
        hbname(id,"lgd_hits",p_lgd,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_lgd_hits:");
      fprintf(stderr," p_lgd pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_lgd_hits - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_time_list(int id, nt_time_list_t *p_tlist)
{
  if (hexist_(&id)) {
    if (p_tlist != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"timelist",p_tlist,NT_FORM_TIME_LIST);
      }
      else {
        hbname(id,"timelist",p_tlist,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_time_list:");
      fprintf(stderr," p_tlist pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_time_list - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_photons(int id, nt_photons_t *p_clust)
{

  if (hexist_(&id)) {
    if (p_clust != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"photons",p_clust,NT_FORM_PHOTONS);
      }
      else {
       hbname(id,"photons",p_clust,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_photons:");
      fprintf(stderr," p_clust pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_photons - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_mesons(int id, nt_mesons_t *p_clust)
{

  if (hexist_(&id)) {
    if (p_clust != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"mesons",p_clust,NT_FORM_MESONS);
      }
      else {
       hbname(id,"mesons",p_clust,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_mesons:");
      fprintf(stderr," p_clust pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_mesons - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_recoil(int id, nt_recoils_t *p_rec)
{

  if (hexist_(&id)) {
    if (p_rec != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"recoils",p_rec,NT_FORM_RECOILS);
      }
      else {
       hbname(id,"recoils",p_rec,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_recoil:");
      fprintf(stderr," p_rec pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_recoil - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_bclusters(int id, nt_bclusters_t *p_bclust)
{

  if (hexist_(&id)) {
    if (p_bclust != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"bcluster",p_bclust,NT_FORM_BCLUSTERS);
      }
      else {
       hbname(id,"bcluster",p_bclust,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_bclusters:");
      fprintf(stderr," p_bclust pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_bclusters - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_tag_coin(int id, nt_tag_coin_t *p_coin)
{

  if (hexist_(&id)) {
    if (p_coin != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"tag_coin",p_coin,NT_FORM_TAG_COIN);
      }
      else {
       hbname(id,"tag_coin",p_coin,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_tag_coin:");
      fprintf(stderr," p_coin pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_tag_coin - ntuple id %d does not exist!\n",
            id);
  }
}

void appendnt_tag_align(int id, nt_tag_align_t *p_talign)
{

  if (hexist_(&id)) {
    if (p_talign != NULL) {
      int entries;
      hnoent(id,&entries);
      if (entries == 0) {
        hbname(id,"tagalign",p_talign,NT_FORM_TAG_ALIGN);
      }
      else {
       hbname(id,"tagalign",p_talign,"$SET");
      }
    }
    else {
      fprintf(stderr,"WARNING: appendnt_tag_align:");
      fprintf(stderr," p_talign pointer was null!\n");
    }
  }
  else {
    fprintf(stderr,"ERROR: appendnt_tag_align - ntuple id %d does not exist!\n",
            id);
  }
}

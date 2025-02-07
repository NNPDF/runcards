{
  if (sd == 1){
    // invariant mass m of the ttx system
    temp_mu_central = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).m();
  }
  else if (sd == 2){
    // transverse mass of the ttx system
    double m  = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).m();
    double pT = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).pT();
    temp_mu_central = sqrt(pow(m, 2) + pow(pT, 2));
  }
  else if (sd == 3){
    // geometric average of transverse masses of and tx - ET
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = sqrt(ET_t * ET_tx);
  }
  else if (sd == 4){
    // sum of transverse masses of t and tx - HT
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = ET_t + ET_tx;
  }
  else if (sd == 5){
    // sum of transverse masses of t and tx and of jets - HTp
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = ET_t + ET_tx;
    for (int i_l = 5; i_l < PARTICLE("all").size(); i_l++){
      temp_mu_central += PARTICLE("all")[i_l].ET;
    }
    // Alternatively:
    // for (int i_l = 3; i_l < PARTICLE("all").size(); i_l++){
  }
  else if (sd == 6){
    // sum of transverse masses of t and tx - HTint
    double xET_t = sqrt(msi->M2_t / 4 + PARTICLE("top")[0].pT2);
    double xET_tx = sqrt(msi->M2_t / 4 + PARTICLE("atop")[0].pT2);
    temp_mu_central = xET_t + xET_tx;
  }
  else if (sd == 7){
    // transverse mass of t - mT_t
    temp_mu_central = PARTICLE("top")[0].ET;
  }
  else if (sd == 8){
    // transverse mass of t - mT_tx
    temp_mu_central = PARTICLE("atop")[0].ET;
  }

  else if (sd == 9){
    // transverse mass of the leading top-jet - mT_t1
    temp_mu_central = PARTICLE("tjet")[0].ET;
  }
  else if (sd == 10){
    // transverse mass of the subleading top-jet - mT_t2
    temp_mu_central = PARTICLE("tjet")[1].ET;
  }
  else if (sd == 11){
    // sum of transverse energy of t and tx at central mtop- HT / 4
    double pT2_t = PARTICLE("top")[0].pT2;
    double pT2_tx = PARTICLE("atop")[0].pT2;
    double ET_t = sqrt(pow(172.5, 2) + pT2_t);
    double ET_tx = sqrt(pow(172.5, 2) + pT2_tx);
    temp_mu_central = 0.25 * (ET_t + ET_tx);
  }
  else if (sd == 70){
    // sum of transverse energies of t and tx varying mtop - HT / 4
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = 0.25 * (ET_t + ET_tx);
  }

  
  /*
  // these were only required because of some former bug in the implementation:
  else if (sd == 11){
    // invariant mass m of the ttx system
    temp_mu_central = (PARTICLE("top")[0].momentum + PARTICLE("atop")[0].momentum).m();
  }
  else if (sd == 14){
    // sum of transverse masses of t and tx - HT
    double ET_t = PARTICLE("top")[0].ET;
    double ET_tx = PARTICLE("atop")[0].ET;
    temp_mu_central = ET_t + ET_tx;
  }
  else if (sd == 17){
    // transverse mass of t - mT_t
    temp_mu_central = PARTICLE("top")[0].ET;
  }
  else if (sd == 18){
    // transverse mass of t - mT_tx
    temp_mu_central = PARTICLE("atop")[0].ET;
  }
  */
  else{
    assert(false && "Scale setting not defined. Reset dynamical scale to different value.");
  }
}

  /*
if (sd == 1){
  fourvector Q;
  for (int i_p = 3; i_p < 5; i_p++){Q = Q + esi->p_parton[i_a][i_p];}
  temp_mu_central = Q.m();
}
  */

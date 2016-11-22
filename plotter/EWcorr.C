double ptWeightQCD(int nGenVbosons, double lheHT, int GenVbosons_pdgId){
  double SF = 1.;
  if (lheHT>100 && nGenVbosons==1){
    if (GenVbosons_pdgId == 23){ // Z
    SF =   ( (lheHT>100 && lheHT<200)*1.588 *(1./1.23 ) + (lheHT>200 && lheHT<400)*1.438 * ( 1./1.23) + (lheHT>400 && lheHT<600)*1.494 * (1./1.23 ) + (lheHT>600)*1.139 * (1./ 1.23) );
    }
    if (abs(GenVbosons_pdgId) == 24){
      SF =   ((lheHT>100 && lheHT<200)* 1.459 * ( 1/ (1.21 ) ) + (lheHT>200 && lheHT<400)* 1.434 * ( 1/ ( 1.21 )) + (lheHT>400 && lheHT<600)*1.532 * (1 / (1.21 )) + (lheHT>600)*1.004 * ( 1 / (1.21) ));
    }
  }
  return SF>0?SF:0;
}

// weights correction for EWK NLO correction
double ptWeightEWK(int nGenVbosons,double GenVbosons_pt,int VtypeSim,int GenVbosons_pdgId){
  double SF = 1.;
  if (nGenVbosons ==1)
    {
      if (VtypeSim == 0 || VtypeSim == 1 || VtypeSim == 4 || VtypeSim == 5)
    {
      if (GenVbosons_pdgId == 23)
        {
          //for Z options
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.1808051+6.04146*(TMath::Power((GenVbosons_pt+759.098),-0.242556));
        }
    }
      else if (VtypeSim == 2 || VtypeSim == 3)
    {
      //for W options
      if (GenVbosons_pdgId == 24 || GenVbosons_pdgId == -24)
        {
          if (GenVbosons_pt > 100. && GenVbosons_pt < 3000) SF = -0.830041+7.93714*(TMath::Power((GenVbosons_pt+877.978),-0.213831));
        }
    }
    }
  return SF>0?SF:0;
}

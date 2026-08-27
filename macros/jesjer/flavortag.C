int JetKinematicCheck::GetJetFlavor(Jet* jet,float jet_radius)
{
  PHHepMCGenEventMap* geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  HepMC::GenEvent* hepmc_event = nullptr;
  if (geneventmap)
    {
      PHHepMCGenEvent* genevt = geneventmap->get(m_embedding_id);
      if (genevt) hepmc_event = genevt->getEvent();
    }

  float jet_eta = jet->get_eta();
  float jet_phi = jet->get_phi();
  int flavor = 0;
  float max_pt = 0;
  for (HepMC::GenEvent::particle_const_iterator p = hepmc_event->particles_begin();
       p != hepmc_event->particles_end(); ++p)
    {
      HepMC::GenParticle *particle = *p;
      if (!particle) continue;
      
      int pid = abs(particle->pdg_id());
      int status = particle->status();
      // Select outgoing partons from hard scattering (status 23) or
      // partons before hadronization (status 21, 22)
      if (status != 23 && status != 21 && status != 22) continue;

      // Only consider quarks (1-6) and gluons (21)
      if (!(pid >= 1 && pid <= 6) && pid != 21) continue;

          HepMC::FourVector momentum = particle->momentum();
          float part_pt = sqrt(momentum.px() * momentum.px() + momentum.py() * momentum.py());
          float part_eta = momentum.eta();
          float part_phi = momentum.phi();

          // Require minimum pT for parton matching
	  if (part_pt < 3.0) continue;

	  // Calculate angular distance between parton and jet
          float dr = CalculateDeltaR(part_eta, part_phi, jet_eta, jet_phi);

	  // Match parton to jet if within jet radius and has highest pT
	  if (dr < jet_radius && part_pt > max_pt)
            {
              max_pt = part_pt;
              flavor = pid;
            }
    }
  return flavor;
}

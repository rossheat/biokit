double dominantPhenotypeProba(
    int countHomozygousDominant, countHeterozygous, countHomozygousRecessive) {
  int totalOrganisms =
      countHomozygousDominant + countHeterozygous + countHomozygousRecessive;

  // Heterozygous * Heterozygous
  double HETHET = ((countHeterozygous / totalOrganisms) *
          ((countHeterozygous - 1) / (totalOrganisms - 1))) *
      0.25;
  // Heterozygous * Homozygous Recessive
  double HETHR = ((countHeterozygous / totalOrganisms) *
      (countHomozygousRecessive / (totalOrganisms - 1)));

  // Homozygous Recessive * Homozygous Recessive
  double HRHR = (countHomozygousRecessive / totalOrganisms) *
      ((countHomozygousRecessive - 1) / (totalOrganisms - 1));

  // 1 - (probability of a recessive phenotype)
  double dominantPhenotypeProbability = 1 - (HETHET + HETHR + HRHR);
  return dominantPhenotypeProbability;
}
